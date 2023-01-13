## ---------------------------
##
## Script Name: Plot LMM Parameters
##
## Purpose of Script: Grpahically display how parameters change with concentration
##
## Author: Zachary Rowson
##
## Date Created: 2022-02-13
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
library(ggplot2)
library(here)

here::i_am("R/plot LMM Parms.R")
model <- lmm1.ARMA.REML

# Extract parameter values and confidence intervals -----------------------

# Extract fixed effects coefficients and CI's
model.int <- intervals(model)
attach(model.int)

fixed.dt <- data.table(fixed,
                       "beta" = c(rep("B_0",10),rep("B_1",10),rep("B_2",10)),
                       "conc" = rep(Padilla_DNT60_Loperamide.L[,unique(conc)],3)
                       )


# Plot fixed effects using ggplot2 ------------------------------------------------------

ggplot() +
  geom_point(fixed.dt[beta=="B_0"], mapping = aes(x=conc,y=`est.`)) +
  geom_errorbar(fixed.dt[beta=="B_0"], mapping = aes(x=conc,ymin=lower,ymax=upper), width = 0.095) +
  labs(title = "lmm1.AR1.REML: Loperamide B_0.D") +
  scale_x_continuous(trans = "log10") +
  theme_bw()

ggplot() +
  geom_point(fixed.dt[beta=="B_1"], mapping = aes(x=conc,y=`est.`)) +
  geom_errorbar(fixed.dt[beta=="B_1"], mapping = aes(x=conc,ymin=lower,ymax=upper), width = 0.095) +
  labs(title = "lmm1.AR1.REML: Loperamide B_1.D") +
  scale_x_continuous(trans = "log10") +
  theme_bw()

ggplot() +
  geom_point(fixed.dt[beta=="B_2"], mapping = aes(x=conc,y=`est.`)) +
  geom_errorbar(fixed.dt[beta=="B_2"], mapping = aes(x=conc,ymin=lower,ymax=upper), width = 0.095) +
  labs(title = "lmm1.AR1.REML: Loperamide B_2.D") +
  scale_x_continuous(trans = "log10") +
  theme_bw()

detach(model.int)


# Calculate and plot var(random effects) ------------------------------------

# Calculate variance in random effects by concentration
z.dt <- as.data.table(model$coefficients$random$fishID, keep.rownames = "fishID")
z.dt2 <- unique(Padilla_DNT60_Loperamide.L[z.dt, on= .(fishID)][, .(fishID,conc,`(Intercept)`)])

var.z.dt <- z.dt2[, .("var" = var(`(Intercept)`)), by = .(conc)]

# Plot variances by concentration
ggplot() +
  geom_point(var.z.dt, mapping = aes(x=conc,y=var)) +
  labs(title = "lmm1.AR1.REML: Loperamide Var(random intercept).L") +
  scale_x_continuous(trans = "log10") +
  theme_bw()
