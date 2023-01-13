## ---------------------------
##
## Script Name: Demonstrate Singular Random Effects Matrix
##
## Purpose of Script: Re-fit random effects with new optimization procedure.
##
## Author: Zachary Rowson
##
## Date Created: 2022-03-04
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


# Format Data -------------------------------------------------------------


load("data/Padilla_DNT60_Loperamide.rda")

# Isolate Loperamide Light data
rmv <- Padilla_DNT60_Loperamide[is.na(y), unique(fishID)]
Padilla_DNT60_Loperamide.L <- Padilla_DNT60_Loperamide[!(fishID%in%rmv) & (t%in%11:30)]

# Specify hierarchical data structure
nest.Lop.L <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_Loperamide.L)

# Set beginning of light to t=0. Set concentration as a factor
attach(nest.Lop.L)
nest.Lop.L$t <- t - min(t)
attach(nest.Lop.)
nest.Lop.L$t2 <- t^2
nest.Lop.L$conc <- as.factor(conc)
detach(nest.Lop.L)
rm(Padilla_DNT60_Loperamide, Padilla_DNT60_Loperamide.L)


# Fit Models and demonstrate singularity ----------------------------------


# Fit base models, random intercept with and w/o variance dependent on concentration
lmm.Lop.L <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                 random = reStruct(~ 0 + conc | fishID, REML=TRUE, pdClass="pdDiag"),
                 data = nest.Lop.L,
                 corr = corAR1(, form = ~ t | fishID), method = "REML")
lmm.Lop.L.noVar <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                       random =  ~ 1 | fishID,
                       data = nest.Lop.L,
                       corr = corAR1(, form = ~ t | fishID), method = "REML")

# Fit random slope term with and w/o variance dependent on concentration
# No randoom intercept term
lmm.Lop.L.rPoly1 <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                        random = reStruct(~ 0 + conc:t | fishID, REML=TRUE, pdClass="pdDiag"),
                        data = nest.Lop.L,
                        corr = corAR1(, form = ~ t | fishID), method = "REML")
lmm.Lop.L.rPoly1.noVar <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                              random = ~ 0 + t | fishID,
                              data = nest.Lop.L,
                              corr = corAR1(, form = ~ t | fishID), method = "REML")

anova(lmm.Lop.L, lmm.Lop.L.noVar, lmm.Lop.L.rPoly1, lmm.Lop.L.rPoly1.noVar)
# Removal of concentration dependent Var(RE) reduces AIC. This may not be the case for all chemicals however.
# Modelling random slope w/o concentration dependent variance is best.

# Include both random intercept and slope.
# Modelling with conc depen. var, we know this doesn't work.
ctrl = lmeControl(maxIter = 200, msMaxIter = 200, minAbsParApVar = 1e-40)
lmm.Lop.L.rIrPoly1 <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                          random = reStruct(~ 0 + conc + conc:t | fishID, REML=TRUE, pdClass="pdDiag"),
                          data = nest.Lop.L,
                          corr = corAR1(, form = ~ t | fishID), method = "REML",
                          control = ctrl)
lmm.Lop.L.rIrPoly1$apVar
lmm.Lop.L.rIrPoly1.noVar <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                                random = ~ 1 + t | fishID,
                                data = nest.Lop.L,
                                corr = corAR1(, form = ~ t | fishID), method = "REML",
                                control = ctrl)
# This results in singular convergence. Potential evidence that random intercept and slope are too highly correlated.

# Fit with random polynomial terms
lmm.Lop.L.rPoly2 <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                        random = reStruct(~ 0 + conc:t + conc:t2 | fishID, REML=TRUE, pdClass="pdDiag"),
                        data = nest.Lop.L,
                        corr = corAR1(, form = ~ t | fishID), method = "REML",
                        control = ctrl) # This results in non-positive definite
lmm.Lop.L.rPoly2.noVar <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                              random = ~ 0 + t + t2 | fishID,
                              data = nest.Lop.L,
                              corr = corAR1(, form = ~ t | fishID), method = "REML",
                              control = ctrl) # This is the best fit
anova(lmm.Lop.L.rPoly1, lmm.Lop.L.rPoly1.noVar, lmm.Lop.L.rPoly2.noVar) # Best fit is random slopes and quads but no Var(conc)

# Experiment with concentration dependent random effects for either slope or quad
lmm.Lop.L.rPoly2.varT <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                            random = reStruct(~ 0 + conc:t + t2 | fishID, REML=TRUE, pdClass="pdDiag"),
                            data = nest.Lop.L,
                            corr = corAR1(, form = ~ t | fishID), method = "REML",
                            control = ctrl) # This results in non-positive definite matrix
lmm.Lop.L.rPoly2.varT2 <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                             random = reStruct(~ 0 + t + conc:t2 | fishID, REML=TRUE, pdClass="pdDiag"),
                             data = nest.Lop.L,
                             corr = corAR1(, form = ~ t | fishID), method = "REML",
                             control = ctrl) # Gives Error "false convergence"
anova(lmm.Lop.L.rPoly2.varT, lmm.Lop.L.rPoly2.varT2)

# Identify the most ideal model that contains a random effect associated with concentration
lmm.Lop.L.onlyRPoly2 <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                            random = reStruct(~ 0 + conc:t2 | fishID, REML=TRUE, pdClass="pdDiag"),
                            data = nest.Lop.L,
                            corr = corAR1(, form = ~ t | fishID), method = "REML",
                            control = ctrl)
anova(lmm.Lop.L, lmm.Lop.L.rPoly1, lmm.Lop.L.onlyRPoly2) # Best is random slope
anova(lmm.Lop.L.rPoly1, lmm.Lop.L.rPoly2.noVar) # Versus the best fit, best fit is significantly better

# TAKE- AWAY
# If we want to model random effect variance as a function of concentration, and use raw polynomial coefficients, we will need to use only a
# random slope effect. It seems that intercept is too highly correlated with other random effects (quad and slope) resulting in singularity.
# Attempting to approximate parameters for two random effects interacting with concentration results in over-parametrization. I believe this
# is due to there only being one observation for each concentration group.
#
# The best model is one where variance of random effects is not a function of concentration, includes random slope and quadratic terms.
# Inclusion of intercept results in singularity.


# Fitting with orthogonal polynomial coefficients --------------------------


ctrl = lmeControl(opt = "optim")
# Try with orthogonal coefficients, variance not modelled as a function of concentration
lmm.Lop.L.rIrPoly1.orthoPoly <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=FALSE),
                                    random = reStruct(~ 0 + conc + conc:poly(t,degree=1,raw=FALSE) | fishID, REML=TRUE, pdClass="pdDiag"),
                                    data = nest.Lop.L,
                                    corr = corAR1(, form = ~ t | fishID), method = "REML",
                                    control = ctrl)
lmm.Lop.L.rIrPoly1.orthoPoly$apVar # Modelling Var(conc) does not work even with orthogonal coeff.

# Try without Var(conc)
lmm.Lop.L.rIrPoly1.noVar.orthoPoly <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=FALSE),
                                          random = ~ 1 + poly(t,degree=1,raw=FALSE) | fishID,
                                          data = nest.Lop.L,
                                          corr = corAR1(, form = ~ t | fishID),
                                          method = "REML",
                                          control = ctrl)
# Fitting works but other models need to be refit prior to comparison

# ReFit base models, random intercept with and w/o variance dependent on concentration
lmm.Lop.L.ortho <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=FALSE),
                       random = reStruct(~ 0 + factor(conc) | fishID, REML=TRUE, pdClass="pdDiag"),
                       corr = corAR1(, form = ~ t | fishID),
                       method = "REML",
                       control = ctrl,
                       data = nest.Lop.L)
lmm.Lop.L.noVar.ortho <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=FALSE),
                       random =  ~ 1 | fishID,
                       data = nest.Lop.L,
                       corr = corAR1(, form = ~ t | fishID),
                       method = "REML")
anova(lmm.Lop.L.rIrPoly1.noVar.orthoPoly, lmm.Lop.L.ortho, lmm.Lop.L.noVar.ortho) # first model is best by far

# ReFit random slope term, no random intercept term
lmm.Lop.L.rPoly1.ortho <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=FALSE),
                        random = reStruct(~ 0 + factor(conc):poly(t,degree=1,raw=FALSE) | fishID, REML=TRUE, pdClass="pdDiag"),
                        data = nest.Lop.L,
                        corr = corAR1(, form = ~ t | fishID),
                        method = "REML",
                        control = ctrl)

lmm.Lop.L.rPoly1.noVar.ortho <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=FALSE),
                              random = ~ 0 + poly(t,degree=1,raw=FALSE) | fishID,
                              data = nest.Lop.L,
                              corr = corAR1(, form = ~ t | fishID),
                              method = "REML",
                              control = ctrl)

anova(lmm.Lop.L.rIrPoly1.noVar.orthoPoly, lmm.Lop.L.rPoly1.ortho, lmm.Lop.L.rPoly1.noVar.ortho) # First model is still best

# Include random quadratic

# Random polynomials with no random intercept term.
lmm.Lop.L.rPoly2.ortho <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=FALSE),
                        random = reStruct(~ 0 + conc:poly(t,degree=2,raw=FALSE) | fishID, REML=TRUE, pdClass="pdDiag"),
                        corr = corAR1(, form = ~ t | fishID),
                        data = nest.Lop.L,
                        method = "REML",
                        control = ctrl)# Non-positive definite

lmm.Lop.L.rPoly2.noVar.ortho <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=FALSE),
                                    random = ~ 0 + poly(t,degree=2,raw=FALSE) | fishID,
                                    data = nest.Lop.L,
                                    corr = corAR1(, form = ~ t | fishID),
                                    method = "REML",
                                    control = ctrl)

# Compare to previous winning model
anova(lmm.Lop.L.rIrPoly1.noVar.orthoPoly, lmm.Lop.L.rPoly2.noVar.ortho) # first model is still best


# Add quad term to first model
lmm.Lop.L.rIrPoly2.noVar.ortho <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=FALSE),
                              random = ~ 1 + poly(t,degree=2,raw=FALSE) | fishID,
                              data = nest.Lop.L,
                              corr = corAR1(, form = ~ t | fishID),
                              method = "REML",
                              control = ctrl)

anova(lmm.Lop.L.rIrPoly1.noVar.orthoPoly, lmm.Lop.L.rIrPoly2.noVar.ortho) # Model with random effects for all polynomial terms is best

# Identify the most ideal model that contains a random effect associated with concentration and orthogonal polynomials
lmm.Lop.L.onlyRPoly2.ortho <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=FALSE),
                            random = reStruct(~ 0 + conc:poly(t,degree=2,raw=FALSE) | fishID, REML=TRUE, pdClass="pdDiag"),
                            data = nest.Lop.L,
                            corr = corAR1(, form = ~ t | fishID), method = "REML",
                            control = ctrl)
anova(lmm.Lop.L.ortho, lmm.Lop.L.rPoly1.ortho, lmm.Lop.L.onlyRPoly2.ortho) # Best is random intercept alone


# TAKE-AWAY
# If we want to model random effect variance as a function of concentration, and use orthogonal polynomial coefficients, we will need to use only a
# random intercept effect. However, since orthogonal coefficients do not increase the number of concentration depenedent parameters we can model,
# using orthogonal coefficients at all is unattractive.
#
# The best model with orthogonal coefficients is one where variance of random effects is not a function of concentration, includes random intercept
# random slope and quadratic terms. Thus the plus of modelling with orthogonal coefficients is that all polynomial terms can be modelled as random effects.

