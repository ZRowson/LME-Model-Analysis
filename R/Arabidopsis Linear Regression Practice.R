#-------------------------------------------
# Walking through a Regression Example
#   with Arabidopsis Data
#   Zachary Rowson (Rowson.Zachary@epa.gov)
#   Created: 02/04/2022
#   Last Edit: 02/04/2022
#-------------------------------------------

library(lme4)
library(nlme)

data("Arabidopsis")
attach(Arabidopsis)

# Overview of the variables
par(mfrow = c(2,4))
barplot(table(reg), ylab = "Frequency", main = "Region")
barplot(table(popu), ylab = "Frequency", main = "Population")
barplot(table(gen), ylab = "Frequency", las = 2, main = "Genotype")
barplot(table(rack), ylab = "Frequency", main = "Rack")
barplot(table(nutrient), ylab = "Frequency", main = "Nutrient")
barplot(table(amd), ylab = "Frequency", main = "AMD")
barplot(table(status), ylab = "Frequency", main = "Status")
hist(total.fruits, col = "grey", main = "Total fruits", xlab = NULL)

# Transform the three factor variables gen, rack and nutrient
Arabidopsis[,c("gen","rack","nutrient")] <- lapply(Arabidopsis[,c("gen","rack","nutrient")], factor)
str(Arabidopsis)
# Re-attach after correction, ignore warnings
attach(Arabidopsis)
# Add 1 to total fruits, otherwise log of 0 will prompt error
total.fruits <- log(1 + total.fruits)

# Overview of the variables
par(mfrow = c(2,4))
barplot(table(reg), ylab = "Frequency", main = "Region")
barplot(table(popu), ylab = "Frequency", main = "Population")
barplot(table(gen), ylab = "Frequency", las = 2, main = "Genotype")
barplot(table(rack), ylab = "Frequency", main = "Rack")
barplot(table(nutrient), ylab = "Frequency", main = "Nutrient")
barplot(table(amd), ylab = "Frequency", main = "AMD")
barplot(table(status), ylab = "Frequency", main = "Status")
hist(total.fruits, col = "grey", main = "Total fruits", xlab = NULL)

# gen x popu table
table(gen, popu)
# Any NAs?
any(is.na(Arabidopsis))


# initial linear regression
LM <- lm(total.fruits ~ rack + nutrient + amd + status)
summary(LM)
par(mfrow = c(2,2))
plot(LM)

# generalized lest squares fitting
GLM <- gls(total.fruits ~ rack + nutrient + amd + status,
           method = "ML")
summary(GLM)

# why fit with ML? REML comparisons are meaningless in LMMs that differ in their fixed effects
