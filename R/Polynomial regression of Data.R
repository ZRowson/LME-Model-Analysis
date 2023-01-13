# apply polynomial regression to row data

library(gabi)
library(data.table)
library(magrittr)

#### trials ###

#### Loperamide ####
par(mfrow = c(1,2))
#### fit to control means ####
load("~/R/gabi/data/DNT60pmr0.rda")
pmr0 <- as.data.table(DNT60pmr0)

# extract data
grp <- gabi::data_egids(pmr0)[cpid == "Loperamide", unique(egid)]
data <- gabi::data_egids(pmr0)[cpid=="Loperamide" | (cpid=="DMSO"&egid==grp)]

# apply polynomial regression to vehicle control data
test.data <- data[wllt=="v"]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L.0 <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D.0 <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Vehicle Control')
predicted.intervals.0 <- predict(fit.L.0,interval='confidence',level=0.99)

t.L.0 <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L.0,predicted.intervals.0[,1],col='green',lwd=3)
lines(t.L.0,predicted.intervals.0[,2],col='black',lwd=1)
lines(t.L.0,predicted.intervals.0[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L.0),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D.0,interval='confidence',level=0.99)

t.D.0 <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D.0,predicted.intervals[,1],col='green',lwd=3)
lines(t.D.0,predicted.intervals[,2],col='black',lwd=1)
lines(t.D.0,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D.0),means[time>30,V1],col='red', pch=19)

#### end ####

#### fit to 0.004 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==0.004]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 0.004 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 0.012 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==0.012]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 0.012 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 0.04 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==0.04]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 0.04 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 0.12 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==0.12]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 0.12 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 0.4 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==0.4]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 0.4 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 1.2 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==1.2]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 1.2 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 4 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==4]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 4 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 12 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==12]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 12 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 40 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==40]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 40 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### 5-Fluorouracil ####

#### fit to control means ####

# extract data
grp <- gabi::data_egids(pmr0)[cpid == "5-Fluorouracil", unique(egid)]
data <- gabi::data_egids(pmr0)[cpid=="5-Fluorouracil" | (cpid=="DMSO"&egid==grp)]

# apply polynomial regression to vehicle control data
test.data <- data[wllt=="v"]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L.0 <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D.0 <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Vehicle Control')
predicted.intervals.0 <- predict(fit.L.0,interval='confidence',level=0.99)

t.L.0 <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L.0,predicted.intervals.0[,1],col='green',lwd=3)
lines(t.L.0,predicted.intervals.0[,2],col='black',lwd=1)
lines(t.L.0,predicted.intervals.0[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L.0),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D.0,interval='confidence',level=0.99)

t.D.0 <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D.0,predicted.intervals[,1],col='green',lwd=3)
lines(t.D.0,predicted.intervals[,2],col='black',lwd=1)
lines(t.D.0,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D.0),means[time>30,V1],col='red', pch=19)

#### end ####

#### fit to 0.12 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==0.12]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 0.004 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 0.4 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==0.4]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 0.012 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 1.2 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==1.2]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 0.04 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 4 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==4]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 0.12 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 12 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==12]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 0.4 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 40 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==40]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 1.2 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####

#### fit to 120 conc means ####
# apply regression to highest concentration tested
test.data <- data[conc==120]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Loperamide at 4 concentration')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

#### end #####



#### 5,5 diphenylhydandoin ####

#### fit to control ####
grp <- gabi::data_egids(pmr0)[cpid == "5,5-diphenylhydandoin", unique(egid)]
data <- gabi::data_egids(pmr0)[cpid=="5,5-diphenylhydandoin" | (cpid=="DMSO"&egid==grp)]

# apply polynomial regression to vehicle control data
test.data <- data[wllt=="v"]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L.0 <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D.0 <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Vehicle Control')
predicted.intervals.0 <- predict(fit.L.0,interval='confidence',level=0.99)

t.L.0 <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L.0,predicted.intervals.0[,1],col='green',lwd=3)
lines(t.L.0,predicted.intervals.0[,2],col='black',lwd=1)
lines(t.L.0,predicted.intervals.0[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L.0),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D.0,interval='confidence',level=0.99)

t.D.0 <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D.0,predicted.intervals[,1],col='green',lwd=3)
lines(t.D.0,predicted.intervals[,2],col='black',lwd=1)
lines(t.D.0,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D.0),means[time>30,V1],col='red', pch=19)
#### end ####

#### fit to 0.12 uM ####
grp <- gabi::data_egids(pmr0)[cpid == "5,5-diphenylhydandoin", unique(egid)]
data <- gabi::data_egids(pmr0)[cpid=="5,5-diphenylhydandoin" | (cpid=="DMSO"&egid==grp)]

# apply polynomial regression to vehicle control data
test.data <- data[conc==0.12]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Vehicle Control')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)
#### end ####
#### fit to 0.40 uM ####
grp <- gabi::data_egids(pmr0)[cpid == "5,5-diphenylhydandoin", unique(egid)]
data <- gabi::data_egids(pmr0)[cpid=="5,5-diphenylhydandoin" | (cpid=="DMSO"&egid==grp)]

# apply polynomial regression to vehicle control data
test.data <- data[conc==0.40]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Vehicle Control')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)
#### end ####


#### fit to 1.20 uM ####
grp <- gabi::data_egids(pmr0)[cpid == "5,5-diphenylhydandoin", unique(egid)]
data <- gabi::data_egids(pmr0)[cpid=="5,5-diphenylhydandoin" | (cpid=="DMSO"&egid==grp)]

# apply polynomial regression to vehicle control data
test.data <- data[conc==1.20]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Vehicle Control')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)
#### end ####
#### fit to 4.00 uM ####
grp <- gabi::data_egids(pmr0)[cpid == "5,5-diphenylhydandoin", unique(egid)]
data <- gabi::data_egids(pmr0)[cpid=="5,5-diphenylhydandoin" | (cpid=="DMSO"&egid==grp)]

# apply polynomial regression to vehicle control data
test.data <- data[conc==4.00]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Vehicle Control')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)
#### end ####
#### fit to 12.00 uM ####
grp <- gabi::data_egids(pmr0)[cpid == "5,5-diphenylhydandoin", unique(egid)]
data <- gabi::data_egids(pmr0)[cpid=="5,5-diphenylhydandoin" | (cpid=="DMSO"&egid==grp)]

# apply polynomial regression to vehicle control data
test.data <- data[conc==12.00]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Vehicle Control')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)
#### end ####
#### fit to 40.00 uM ####
grp <- gabi::data_egids(pmr0)[cpid == "5,5-diphenylhydandoin", unique(egid)]
data <- gabi::data_egids(pmr0)[cpid=="5,5-diphenylhydandoin" | (cpid=="DMSO"&egid==grp)]

# apply polynomial regression to vehicle control data
test.data <- data[conc==40.00]

## melt data
to.mdl <- data.table::melt(test.data, id.vars=names(test.data)[1:9], measure.vars=names(test.data)[10:59],
                           value.name="speed", variable.name="time")
to.mdl[, time:=gsub("vt","",time)]
to.mdl[, time:=as.numeric(time)]

fit.L <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = to.mdl[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(to.mdl[,time],to.mdl[,speed],col='deepskyblue4',xlab='q',main='Vehicle Control')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- to.mdl[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## add mean values
means <- to.mdl[, mean(speed,na.rm=TRUE), by=time]
points(unique(t.L),means[time>10&time<31,V1],col='red', pch=19)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- to.mdl[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## add mean values
points(unique(t.D),means[time>30,V1],col='red', pch=19)
#### end ####


#### Test ####
#### fit regression to fish individually ####

# fit regression to each fish and show distribution of regression parameters

# set up data
to.mdl[, fish := paste0(apid,rowi,coli)]
fish.ids <- to.mdl[, unique(fish)]

#### trial run ####
fish.id <- fish.ids[1]
fish.data <- to.mdl[fish==fish.id]

fit.L <- lm(speed ~ poly(time,2), data = fish.data[!is.na(speed) & time%in%c(11:30)])
fit.D <- lm(speed ~ poly(time,2), data = fish.data[!is.na(speed) & time%in%c(31:50)])

## plot fit in Light
plot(fish.data[,time],fish.data[,speed],col='deepskyblue4',xlab='t',main='One Control Fish')
predicted.intervals <- predict(fit.L,interval='confidence',level=0.99)

t.L <- fish.data[!is.na(speed) & time%in%c(11:30), time]
lines(t.L,predicted.intervals[,1],col='green',lwd=3)
lines(t.L,predicted.intervals[,2],col='black',lwd=1)
lines(t.L,predicted.intervals[,3],col='black',lwd=1)

## plot fit in Dark
predicted.intervals <- predict(fit.D,interval='confidence',level=0.99)

t.D <- fish.data[!is.na(speed) & time>30, time]
lines(t.D,predicted.intervals[,1],col='green',lwd=3)
lines(t.D,predicted.intervals[,2],col='black',lwd=1)
lines(t.D,predicted.intervals[,3],col='black',lwd=1)

## extract parameter values
parms <- fit.D$coefficients

#### end ####


#### apply to all fish ####
fish.ids <- to.mdl[!is.na(speed), unique(fish)]
parms <- sapply(fish.ids, function(id) {

  # identify data
  data <- to.mdl[fish==id]

  # fit regression
  fit.L <- lm(speed ~ poly(time,2), data = data[!is.na(speed) & time%in%c(11:30)])
  fit.D <- lm(speed ~ poly(time,2), data = data[!is.na(speed) & time%in%c(31:50)])

  ## extract parameter values
  list(fit.L$coefficients, fit.D$coefficients)
})

parms <- lapply(apply(parms, 2, rbind), unlist)
parms <- do.call('rbind', parms)
colnames(parms) <- c("B_0.L","B_1.L","B_2.L",
                     "B_0.D","B_1.D","B_2.D")

# compare average coefficients to acceleration values

## calculate averages
mean.parms <- apply(parms, 2, mean, na.rm = TRUE)

# ## plot
# func.L <- function(x) {mean.parms[1] + scale*mean.parms[2]*x + scale*mean.parms[3]*(x^2)}
# func.D <- function(x) {mean.parms[4] + mean.parms[5]*x + mean.parms[6]*(x^2)}
#
# x.L <- seq(0,19,by=0.01)
# x.D <- seq(31,50,by=0.01)
# lines(seq(11,30,by=0.01), scale*func.L(x.L),col='red',lwd=3)
# lines(seq(31,50,by=0.01),predicted.intervals[,2],col='black',lwd=1)
