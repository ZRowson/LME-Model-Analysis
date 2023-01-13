## ---------------------------
##
## Script Name: Implement Woody's suggestions
##
## Purpose of Script:
##
## Author: Zachary Rowson
##
## Date Created: 2022-03-15
##
## Email: Rowson.Zachary@epa.gov
##
## Working Directory:
##
## ---------------------------
##
## Notes:
## Woody and I met 03/15/2022 initially to discuss optimizing the zebrafish protocol.
##      *Woody states he managed to get the model to fit by changing the default apVar structure in the lme call
##        * apparently the default is not pdsymm but is a log-cholinksy something or other
##        * somewhere along the line this matrix likely explodes
##        * Woody is also trying to fit models to the log transformed data but finds that some random effects become obsolete
##
## ---------------------------


library(nlme)
library(data.table)
library(ggplot2)


load("data/Padilla_DNT60_Loperamide.rda")

data <- copy(Padilla_DNT60_Loperamide)
data.L <- Padilla_DNT60_Loperamide[t%in%11:30]
data.D <- Padilla_DNT60_Loperamide[t%in%31:50]
remove <- Padilla_DNT60_Loperamide[is.na(y), unique(fishID)]

nest.L <- groupedData(y ~ conc|fishID, data = data.L[!(fishID%in%remove)], order.groups = FALSE)
nest.D <- groupedData(y ~ conc|fishID, data = data.D[!(fishID%in%remove)], order.groups = FALSE)

# Set beginning of light to t=0. Set concentration as a factor
nest.L$t <- nest.L$t - min(nest.L$t)
nest.L$conc <- as.factor(nest.L$conc)
nest.D$t <- nest.D$t - min(nest.D$t)
nest.D$conc <- as.factor(nest.D$conc)


# Fit pdSymm model to Light --------------------------------------------------------


# Light
ctrl <- lmeControl(opt = "optim", niterEM = 200, msVerbose = TRUE, msTol = 1e-20, msMaxIter = 200)
lmm1 <- lme(y ~  conc*poly(t,degree=2,raw=TRUE),
            correlation = corAR1(, form = ~ t | fishID),
            random = reStruct(~1|fishID,REML=FALSE,pdClass="pdSymm"),
            method = "ML",
            control = ctrl,
            data = nest.L)
lmm2 <- lme(y ~ conc*poly(t,degree=2,raw=TRUE),
            correlation = corAR1(, form = ~ t | fishID),
            random = reStruct(~poly(t,degree=1,raw=TRUE)|fishID,REML=FALSE,pdClass="pdSymm"),
            method = "ML",
            control = ctrl,
            data = nest.L)
lmm3 <- lme(y ~ conc*poly(t,degree=2,raw=TRUE),
            correlation = corAR1(, form = ~ t | fishID),
            random = reStruct(~poly(t,degree=2,raw=TRUE)|fishID,REML=FALSE,pdClass="pdSymm"),
            method = "ML",
            control = ctrl,
            data = nest.L)
lmm1$apVar
lmm2$apVar
lmm3$apVar
anova(lmm1,lmm2, lmm3) # lmm3 is best


# Power transform data and apply LMEs -------------------------------------


# Save log10 data data.table
data.L[,logy:=log10(y+1)]
data.D[,logy:=log10(y+1)]

# Set up hierarchical datasets for light and dark.
nest.log.L <- groupedData(logy ~ conc|fishID, data = data.L[!(fishID%in%remove)], order.groups = FALSE)
nest.log.D <- groupedData(logy ~ conc|fishID, data = data.D[!(fishID%in%remove)], order.groups = FALSE)

# Light
ctrl <- lmeControl(opt = "optim", niterEM = 200, msVerbose = TRUE, msTol = 1e-10, msMaxIter = 400, allow.n.lt.q = TRUE)
lmm1.log <- lme(logy ~  conc*poly(t,degree=2,raw=TRUE),
            correlation = corAR1(, form = ~ t | fishID),
            random = reStruct(~1|fishID,REML=FALSE,pdClass="pdSymm"),
            method = "ML",
            control = ctrl,
            data = nest.log.L)
lmm2.log <- lme(logy ~ conc*poly(t,degree=2,raw=TRUE),
            correlation = corAR1(, form = ~ t | fishID),
            random = reStruct(~poly(t,degree=1,raw=TRUE)|fishID,REML=FALSE,pdClass="pdSymm"),
            method = "ML",
            control = ctrl,
            data = nest.log.L)
lmm3.log <- lme(logy ~ conc*poly(t,degree=2,raw=TRUE),
            correlation = corAR1(, form = ~ t | fishID),
            random = reStruct(~poly(t,degree=2,raw=TRUE)|fishID,REML=FALSE,pdClass="pdSymm"),
            method = "ML",
            control = ctrl,
            data = nest.log.L)
lmm1$apVar
lmm2$apVar
lmm3$apVar
anova(lmm1.log, lmm3.log)

# Play with models in dark ------------------------------------------------



ctrl <- lmeControl(opt = "optim", niterEM = 200, msVerbose = TRUE, msTol = 1e-20, msMaxIter = 200)
lmm1D <- lme(y ~  0 + conc*poly(t,degree=2,raw=TRUE),
            correlation = corARMA(,form=~t|fishID,p=1,q=1),
            random = reStruct(~1|fishID,REML=FALSE,pdClass="pdSymm"),
            method = "ML",
            data = nest.D) # optim optimizer fails when using ARMA "non-invertible coefficient matrix", singularity?
lmm2D <- lme(y ~ conc*poly(t,degree=2,raw=TRUE),
            correlation = corARMA(,form=~1|fishID,p=1,q=1),
            random = reStruct(~poly(t,degree=1,raw=TRUE)|fishID,REML=FALSE,pdClass="pdSymm"),
            method = "ML",
            control = ctrl,
            data = nest.L)
lmm3D <- lme(y ~ conc*poly(t,degree=2,raw=TRUE),
            correlation = corARMA(,form=~t|fishID,p=1,q=1),
            random = reStruct(~poly(t,degree=2,raw=TRUE)|fishID,REML=FALSE,pdClass="pdSymm"),
            method = "ML",
            control = ctrl,
            data = nest.L)
lmm1$apVar
lmm2$apVar
lmm3$apVar
anova(lmm1,lmm2, lmm3) # lmm3 is best
