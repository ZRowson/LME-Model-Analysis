## ---------------------------
##
## Script Name: Re-fit Random Effects
##
## Purpose of Script: Re-fit random effects with new optimization procedure.
##
## Author: Zachary Rowson
##
## Date Created: 2022-02-21
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
library(gabi)


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
rm(DNT60pmr0_long, grp)

# Specify hierarchical data structure
nest.55.L <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_55diphenyl.L[!(fishID%in%rmv)])
attach(nest.55.L)
nest.55.L$t <- t - min(t)
nest.55.L$conc <- as.factor(conc)
detach(nest.55.L)


rm(Padilla_DNT60_55diphenyl.L, Padilla_DNT60_Loperamide, Padilla_DNT60_Loperamide.L)


# Experiment with Loperamide in Light -------------------------------------


# Fit base model, random intercept with variance changing with concentration
lmm.Lop.L <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                 random = reStruct(~ 0 + factor(conc) | fishID, REML=TRUE, pdClass="pdDiag"),
                 data = nest.Lop.L,
                 corr = corAR1(, form = ~ t | fishID), method = "REML")

ggplot(nest.Lop.L, aes(t,y, color = conc)) +
  geom_point()+
  geom_line(aes(group=fishID))+
  facet_grid(.~conc)

# Fit random slope term with and without Random Effect variance dependent on concentration
lmm.Lop.L.rPoly1 <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                 random = reStruct(~ 0 + factor(conc):poly(t,degree=1,raw=TRUE) | fishID, REML=TRUE, pdClass="pdDiag"),
                 data = nest.Lop.L,
                 corr = corAR1(, form = ~ t | fishID), method = "REML")

lmm.Lop.L.rPoly1.noVar <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                        random = ~ 0 + poly(t,degree=1,raw=TRUE) | fishID,
                        data = nest.Lop.L,
                        corr = corAR1(, form = ~ t | fishID), method = "REML")

anova(lmm.Lop.L, lmm.Lop.L.rPoly1, lmm.Lop.L.rPoly1.noVar) # Removal of concentration dependent Var(RE) reduces AIC

# Include both random intercept and slope
ctrl = lmeControl(maxIter = 200, msMaxIter = 200)
lmm.Lop.L.rIrPoly1.noVar <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                              random = ~ 1 + poly(t,degree=1,raw=TRUE) | fishID,
                              data = nest.Lop.L,
                              corr = corAR1(, form = ~ t | fishID), method = "REML",
                              control = ctrl) # This results in singular convergence

# Try with orthogonal coefficients
ctrl = lmeControl(opt = "optim")
lmm.Lop.L.rIrPoly1.noVar.noRaw <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=FALSE),
                                random = ~ 1 + poly(t,degree=1,raw=FALSE) | fishID,
                                data = nest.Lop.L,
                                corr = corAR1(, form = ~ t | fishID), method = "REML",
                                control = ctrl) # Fitting works but other models need to be refit prior to comparison




# Fit random quadratic term and compare AIC
lmm.Lop.L.rPoly2 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest.Lop.L, method = "REML",
                          corr = corAR1(,form=~t|fishID),
                          random = reStruct(~ 0 + factor(conc):poly(t,degree=2,raw=TRUE) | fishID, REML=TRUE, pdClass="pdDiag"),
                          control = ctrl)

lmm.Lop.L.rPoly2.noVar <- lme(y ~ 0 + conc + conc:poly(t,degree=2, raw=TRUE),
                              random = ~ 0 + poly(t,degree=2,raw=TRUE) | fishID,
                              data = nest.Lop.L,
                              corr = corAR1(, form = ~ t | fishID), method = "REML")

anova(lmm.Lop.L, lmm.Lop.L.rPoly1, lmm.Lop.L.rPoly1.noVar, lmm.Lop.L.rPoly2, lmm.Lop.L.rPoly2.noVar)



# Final model will be a LME with random slope term with Var dependent on concentration
lmm.final.Lop.L <- lmm1.AR1.REML.rPoly1.var.Lop.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE),
                                                         random = reStruct(~ 0 + factor(conc):poly(t,degree=1,raw=TRUE) | fishID, REML = TRUE, pdClass="pdDiag"),
                                                         data = nest.Lop.L,
                                                         corr = corAR1(, form = ~ t | fishID), method = "REML")
intervals(lmm.final.Lop.L)
