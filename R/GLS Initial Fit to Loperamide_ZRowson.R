#----------------------------------------------------
# Apply Generalized Least Squares to a Trial Dataset
#   Zachary Rowson (Rowson.Zachary@epa.gov)
#   Created: 02/03/2022
#   Last Edit: 02/03/2022
#----------------------------------------------------

library(nlme)
library(data.table)

# load data
load("data/Padilla_DNT60_Loperamide.rda")

# remove incomplete data
rmv <- Padilla_DNT60_Loperamide[is.na(y), unique(fishID)]
to.fit <- Padilla_DNT60_Loperamide[!(fishID%in%rmv)]


#### fit to Light data #####

# provide hierarchical structure
nestinginfo.L <- groupedData(y ~ conc|fishID, data =  to.fit[t%in%11:30])

# fit gls with AR(1) correlation structure and polynomial 2 term
fit.ar1polytime.L <- gls(y ~ factor(conc) * poly(t,degree=2), data = nestinginfo.L, corr = corAR1(, form=~1|fishID))

#### end ####

#### fit to Light data ####

# provide hierarchical structure
nestinginfo.D <- groupedData(y ~ conc|fishID, data =  to.fit[t%in%31:50])

# fit gls with AR(1) correlation structure and polynomial 2 term
fit.ar1polytime.D <- gls(y ~ factor(conc) * poly(t,degree=2), data = nestinginfo.D, corr = corAR1(, form=~1|fishID))

#### end ####

summary(fit.ar1polytime.L)
summary(fit.ar1polytime.D)
