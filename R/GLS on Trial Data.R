#---------------------------------
# Apply Generalized Least Squares
# to a trial dataset
#   Zachary Rowson
#   Rowson.Zachary@epa.gov
#---------------------------------

library(nlme)
library(data.table)
library(magrittr)

phlebitisdata <- read.table("phlebitis.csv", header=T, sep=",")
attach(phlebitisdata)

# create hierarchical structure
nestinginfo <- groupedData(Y ~ Treatment | Animal, data = phlebitisdata)

# fit regression with poly 2 structure and AR1 correlations of repeated measures
# lag is calculated automatically, phi = 0.2403242 (Is this the same as rho?)
fit.ar1polytime <- gls(Y ~ factor(Treatment)*poly(Time, degree = 2),
                       data=nestinginfo, corr=corAR1(,form= ~ 1 | Animal))

# try on a sample chemical

# isolate data
chm <- "Loperamide"
grp <- gabi::data_egids(pmr0)[cpid==chm, unique(egid)]
data <- gabi::data_egids(pmr0)[cpid==chm | (wllt=="v"&egid==grp)]

# elongate data
data[, fishID := paste(apid,rowi,coli,sep=".")]
data_long <- data.table::melt(data, id.vars = c(names(data)[1:9],"fishID"), measure.vars = names(data)[10:59],
                              value.name = "y", variable.name = "t")
data_long[, t := as.numeric(gsub("vt","",t))]

# remove blank entries
rmv <- data_long[is.na(y), unique(fishID)]
to.fit <- data_long[!(fishID%in%rmv) & (t%in%11:30)]
# provide hierarchical structure
nestinginfo <- groupedData(y ~ conc|fishID, data =  to.fit)
fit.ar1polytime.L <- gls(y ~ factor(conc) * poly(t,degree=2), data = nestinginfo, corr = corAR1(, form=~1|fishID))

