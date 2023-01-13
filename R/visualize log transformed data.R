## ---------------------------
##
## Script Name: Visualize log10 data for various chemicals
##
## Purpose of Script: Visualize transformed data
##
## Author: Zachary Rowson
##
## Date Created: 2022-03-17
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


library(ggplot2)
library(data.table)


# Load data and transform


load("data/Padilla_DNT60_pmr0_long.rda")
DNT60pmr0_long[, logy := log10(y + 1)]
DNT60pmr0_long[, sqrty := sqrt(y)]

# Display log transformed data over all chemicals
ggplot(DNT60pmr0_long, aes(x=factor(t),y=logy)) +
  geom_boxplot()
# Display sqrt transformed data over all chemicals
ggplot(DNT60pmr0_long, aes(x=factor(t),y=sqrty)) +
  geom_boxplot()
# Visually sqrt transformed data appears to have lower degree of heteroscedasticity in the Dark

# Perform very poor visual assessment of the heteroscedasticity of the data
powerSD <- DNT60pmr0_long[, .("SD" = stats::sd(y,na.rm=TRUE),
                         "logSD" = stats::sd(logy,na.rm=TRUE),
                         "sqrtSD" = stats::sd(sqrty,na.rm=TRUE)), by = t]
powerCV <- DNT60pmr0_long[, .("SD" = stats::sd(y,na.rm=TRUE)/mean(y,na.rm=TRUE),
                              "logSD" = stats::sd(logy,na.rm=TRUE)/mean(logy,na.rm=TRUE),
                              "sqrtSD" = stats::sd(sqrty,na.rm=TRUE)/mean(sqrty,na.rm=TRUE)), by = t]
powerSDlong <- melt(powerSD, id.vars = "t")
powerCVlong <- melt(powerCV, id.vars = "t")
powerSDlong[t%in%11:30, phase :="Light"][t%in%31:50, phase := "Dark"]
powerCVlong[t%in%11:30, phase :="Light"][t%in%31:50, phase := "Dark"]

# Plot histograms of standard deviation data
ggplot(powerSDlong[variable=="SD"], aes(x=value)) +
  geom_histogram() +
  facet_wrap(vars(phase))
ggplot(powerSDlong[variable=="logSD"], aes(x=value)) +
  geom_histogram() +
  facet_wrap(vars(phase))
ggplot(powerSDlong[variable=="sqrtSD"], aes(x=value)) +
  geom_histogram() +
  facet_wrap(vars(phase))

# From these graphics it apears that log10 transform actually reduces heteroscedasticity the most in the Light
# sqrt is better for Dark as expected. However, the scale of these values is much smaller for log and the magnitude of the
# difference in variances will be smaller. Maybe this is all the model cares about

# Maybe log10 transformed data appears the nicest but this is somewhat unclear

# Display histograms of lowest activity time periods
hist(DNT60pmr0_long[t==11,y])
hist(DNT60pmr0_long[t==11,logy])
hist(DNT60pmr0_long[t==11,sqrty])

hist(DNT60pmr0_long[t==31,y])
hist(DNT60pmr0_long[t==31,logy])
hist(DNT60pmr0_long[t==31,sqrty])

# From these two histograms it appears that sqrt transform is best
# I think the deciding factor is going to be an evaluation of the residuals of the data as well as some kind
# of R^2 metric or measure of the variability explained by the model.

ggplot(DNT60pmr0_long[t==11&y!=0], aes(x=sqrty)) +
  geom_histogram(position="identity", alpha = 0.4, binwidth = 0.01) +
  labs(title = "Distribution of Square Root of Total Movement from 20-22 minutes",
       x = "Square Root of Total Distance (sqrt(cm))",
       fill = "time period")
