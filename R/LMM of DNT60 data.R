#-------------------------------------------
# LMM on DNT60 Data
#   Zachary Rowson (Rowson.Zachary@epa.gov)
#   Created: 02/07/2022
#   Last Edit: 02/08/2022
#-------------------------------------------

library(nlme)
library(lme4)
library(data.table)

# load data
load("data/Padilla_DNT60_Loperamide.rda")

# Isolate Dark data first
rmv <- Padilla_DNT60_Loperamide[is.na(y), unique(fishID)]
Padilla_DNT60_Loperamide.D <- Padilla_DNT60_Loperamide[!(fishID%in%rmv) & (t%in%31:50)]
Padilla_DNT60_Loperamide.D[, conc := as.character(conc)]

# Specify hierarchical structure of repeated measures
nestData.D <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_Loperamide.D)

# Shift time variable to zero
nestData.D$t <- nestData.D$t - min(nestData.D$t)

# Fit data with GLS with various fixed effects
lm1 <- gls(y ~ factor(as.character(conc)) * t, data = nestData.D, method = "ML")
lm2.poly2 <- gls(y ~ factor(conc) * poly(t,degree=2,raw=TRUE), data = nestData.D, method = "ML")
lm2.poly3 <- gls(y ~ factor(conc) * poly(t,degree=3,raw=TRUE), data = nestData.D, method = "ML")

anova(lm1, lm2.poly2, lm2.poly3)

# Possibility that cubic term is unnecessary. Try to include plate number
lm3.poly2 <- update(lm2.poly2, .~. + apid)
lm4.poly2 <- update(lm2.poly2, .~. + apid + rowi + coli)
lm3.poly3 <- update(lm2.poly3, .~. + apid)
lm4.poly3 <- update(lm2.poly3, .~. + apid + rowi + coli)

anova(lm2.poly2, lm2.poly3, lm3.poly2, lm3.poly3, lm4.poly2, lm4.poly3) # it seems that necessity of poly3 term is still questionable
# lm4.poly2 seems to be most ideal model
# problem with including apid rowi and coli is that chemical might be nested in these things....

# add a correlation structure
lm6 <- update(lm5, .~., corr = corAR1(, form = ~1|fishID))
lm6.2 <- update(lm5.2, .~., corr = corAR1(, form = ~1|fishID))
lm7 <- update(lm5, .~., corr = corARMA(,form = ~1|fishID,p=1,q=1))
lm7.2 <- update(lm5.2, .~., corr = corARMA(,form = ~1|fishID,p=1,q=1))
anova(lm5, lm5.2, lm6, lm6.2, lm7, lm7.2)

# ARMA(p=1, q=1) seems to be best

# Add random effects but with no covariance structure
lmm1 <- lme(y ~ apid + rowi + coli + apid + factor(conc) * poly(t, degree=2),
            random = ~1|fishID, method = "ML",
            data = nestData.D)
plot(ranef(lmm1))
lmm2 <- lme(y ~ apid + rowi + coli + apid + factor(conc) * poly(t, degree=2),
            random = ~1|apid/fishID, method = "ML",
            data = nestData.D)
plot(ranef(lmm2)) # random effects look much more normal

anova(lmm1, lmm2) # AIC doesn't change

# Add correlation structure. ARMA wins again
lmm3 <- update(lmm1, corr = corARMA(,form = ~1|fishID,p=1,q=1))
lmm4 <- update(lmm1, corr = corAR1(,form = ~1|fishID))
anova(lmm3, lmm4)

lme(y ~ rowi + coli + apid + factor(conc) * poly(t, degree=2),
    random = ~1|fishID, method = "ML", corr = corARMA(,form=~1|fishID,p=1,q=1),
    data = nestData.D)

nestData.D2 <- groupedData(y ~ conc|apid/fishID, data =  Padilla_DNT60_Loperamide.D)
