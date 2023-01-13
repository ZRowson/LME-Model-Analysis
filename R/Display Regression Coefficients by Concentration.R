#--------------------------------------------------
# Display Regression Coefficients by Concentration
#   Zachary Rowson (Rowson.Zachary@epa.gov)
#   Created: 01/27/2022
#   Last Edit: 01/31/2022
#--------------------------------------------------

library(data.table)
library(ggplot2)

# First attempt look at 5,5-diphenylhydandoin. Gradual concentration dependent
# change in slope and curvature in Dark seems obvious in this example
chm <- "5,5-diphenylhydandoin"
gabi::plot_tSeries(pmr0, chemical = chm, prsp = "SS")

#### apply poly2 regression to each concentration group ####

# load data
load("~/R/gabi/data/DNT60pmr0.rda")
pmr0 <- as.data.table(DNT60pmr0)

# identify concentration groups
grp <- gabi::data_egids(pmr0)[cpid==chm, unique(egid)]
concs <- gabi::data_egids(pmr0)[cpid==chm | (wllt=="v"&egid==grp), unique(conc)]

# extract relevant chemical exposure data, elongate, and format
data <- gabi::data_egids(pmr0)[cpid==chm | (wllt=="v"&egid==grp)]
data_long <- data.table::melt(data, id.vars=names(data)[1:9], measure.vars=names(data)[10:59],
                 value.name="speed", variable.name="time")
data_long[, time := as.numeric(gsub("vt","",time))]

## examine sample sizes by concentration group
data[,.N, by=conc]
data[wllq==1, .N, by=conc]

## identify fish data that are incomplete
rmv <- data[rowSums(is.na(data)) >0, names(data)[1:9], with=F]


# fit poly2 to data
fits <- lapply(concs, function(cnc) {

  to.fit <- data_long[!rmv, on=names(rmv)][conc==cnc]

  fit.L <- lm(speed ~ poly(time,2), data = to.fit[time%in%c(11:30)])
  fit.D <- lm(speed ~ poly(time,2), data = to.fit[time%in%c(31:50)])

  list("fit.L" = fit.L, "fit.D" = fit.D)
})
names(fits) <- concs

# parms <- lapply(apply(parms, 2, rbind), unlist)
# parms <- do.call('rbind', parms)
# colnames(parms) <- c("B_0.L","B_1.L","B_2.L",
#                      "B_0.D","B_1.D","B_2.D")

#### end ####

#### plot regression fits ####

# extract coefficient data for plotting. Lots of steps being squished together here, sorry
coeff <- lapply(fits, lapply, `[[`, "coefficients")
coeff.dt <- as.data.table(do.call('rbind', lapply(coeff, unlist)), keep.rownames = "conc")
data.table::setnames(coeff.dt, names(coeff.dt)[-1], c("B_0.L","B_1.L","B_2.L","B_0.D","B_1.D","B_2.D"))

# get confidence intervals for plotting and edit
CIs <- lapply(fits, lapply, confint)
CIs.dt <- as.data.table(do.call('rbind', lapply(CIs,do.call,what='rbind')), keep.rownames = "coeff")
CIs.dt[coeff=="(Intercept)", coeff := "B_0"]
CIs.dt[coeff=="poly(time, 2)1", coeff := "B_1"]
CIs.dt[coeff=="poly(time, 2)2", coeff := "B_2"]
CIs.dt[, conc := as.numeric(rep(names(fits),each=6))]
CIs.dt[, phase := rep(c("Light","Dark"), each=3, times=length(concs))]
data.table::setnames(CIs.dt, names(CIs.dt)[c(2,3)], c("min","max"))

# plot

## Dark
ggplot() +
  geom_point(coeff.dt, mapping = aes(x=conc,y=B_1.D)) +
  geom_errorbar(CIs.dt[coeff=="B_1" & phase=="Dark"], mapping = aes(x=conc,ymin=min,ymax=max), width = 0.095) +
  labs(title = "5,5-diphenylhydandoin: B_1.D") +
  scale_x_continuous(trans = "log10") +
  theme_bw()

ggplot() +
  geom_point(coeff.dt, mapping = aes(x=conc,y=B_2.D)) +
  geom_errorbar(CIs.dt[coeff=="B_2" & phase=="Dark"], mapping = aes(x=conc,ymin=min,ymax=max), width = 0.095) +
  labs(title = "5,5-diphenylhydandoin: B_2.D") +
  scale_x_continuous(trans = "log10") +
  theme_bw()

## Light
ggplot() +
  geom_point(coeff.dt, mapping = aes(x=conc,y=B_1.L)) +
  geom_errorbar(CIs.dt[coeff=="B_1" & phase=="Light"], mapping = aes(x=conc,ymin=min,ymax=max), width = 0.095) +
  labs(title = "5,5-diphenylhydandoin: B_1.L") +
  scale_x_continuous(trans = "log10") +
  theme_bw()

ggplot() +
  geom_point(coeff.dt, mapping = aes(x=conc,y=B_2.L)) +
  geom_errorbar(CIs.dt[coeff=="B_2" & phase=="Light"], mapping = aes(x=conc,ymin=min,ymax=max), width = 0.095) +
  labs(title = "5,5-diphenylhydandoin: B_2.L") +
  scale_x_continuous(trans = "log10") +
  theme_bw()
#### end ####

