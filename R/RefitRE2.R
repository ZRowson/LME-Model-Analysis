## ---------------------------
##
## Script Name: Re-Fit Random Effects 2
##
## Purpose of Script: Fit random effects following Woody's code.
##
## Author: Zachary Rowson
##
## Date Created: 2022-03-11
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
library(nlme)


# Load data ---------------------------------------------------------------


load("data/Padilla_DNT60_Loperamide.rda")

# Isolate Loperamide Light data
rmv <- Padilla_DNT60_Loperamide[is.na(y), unique(fishID)]
Padilla_DNT60_Loperamide.L <- Padilla_DNT60_Loperamide[!(fishID%in%rmv) & (t%in%11:30)]

# Specify hierarchical data structure
nest.Lop.L <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_Loperamide.L)
attach(nest.Lop.L)
nest.Lop.L$t <- t - min(t)
detach(nest.Lop.L)
attach(nest.Lop.L)
nest.Lop.L$t2 <- t^2
nest.Lop.L$conc <- as.factor(conc)
detach(nest.Lop.L)

rm(Padilla_DNT60_Loperamide, Padilla_DNT60_Loperamide.L)


# Create function to normalize AIC ----------------------------------------


# Function orders models by AIC after subtracting minimum AIC observed in set of models being compared
normAIC <- function(zz) {
  zz$AIC <- zz$AIC - min(zzz$AIC)
  indx <- order(zz$AIC)
  zz[indx,]
}


# Start fitting models ----------------------------------------------------
ctrl = lmeControl(maxIter = 200, msMaxIter = 200, minAbsParApVar = 1e-40, opt = "optim") # I believe that first three arguments become obsolete once we specify opt = "optim"

# Best fitting models fit by Woody
lmm.Lop.L.rIrPoly2.noVar.noVar <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
                                      random = reStruct(~ 1 + t + t2 | fishID, REML=FALSE, pdClass="pdDiag"),
                                      corr = corAR1(, form = ~ t | fishID),
                                      method = "ML",
                                      control = ctrl,
                                      data = nest.Lop.L)

lmm.Lop.L.rIrPoly2.noVar.noVar.full  <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
                                            random = reStruct(~ poly(t,degree=2,raw=TRUE) | fishID, REML=FALSE),
                                            corr = corAR1(, form = ~ t | fishID),
                                            method = "ML",
                                            control = ctrl,
                                            data = nest.Lop.L)

t <- seq(from=11, to=20, by = 1)
y <-
