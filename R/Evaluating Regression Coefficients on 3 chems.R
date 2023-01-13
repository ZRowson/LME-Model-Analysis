#-------------------------------------------
# Evaluating Regression Coefficients
# Fit to Three Chemicals
#   Zachary Rowson (Rowson.Zachary@epa.gov)
#   Created: 01/31/2022
#   Last Edit: 01/31/2022
#-------------------------------------------

library(data.table)

# load pmr0 data
load("~/R/gabi/data/DNT60pmr0.rda")
pmr0 <- as.data.table(DNT60pmr0)

# identify chemicals
chms <- c("Loperamide", "5-Fluorouracil", "5,5-diphenylhydandoin")

#### visualize behavior data by concentration for each chemical ####
gabi::plot_tSeries(pmr0, chemical = chms[1], prsp = "SS")
gabi::plot_tSeries(pmr0, chemical = chms[2], prsp = "SS")
gabi::plot_tSeries(pmr0, chemical = chms[3], prsp = "SS")
#### end ####

#### produce poly2 regression fits for each chemical at each concentration ####
fits <- lapply(chms, function(chm) {

  # find concentrations associated with chemical and corresponding vehicle control
  pmr0_w_egid <- gabi::data_egids(pmr0)
  grp <- pmr0_w_egid[cpid==chm, unique(egid)]
  concs <- pmr0_w_egid[cpid==chm | (wllt=="v"&egid==grp), unique(conc)]

  # extract relevant behavioral data, elongate, and format
  data <- pmr0_w_egid[cpid==chm | (wllt=="v"&egid==grp)]
  data_long <- data.table::melt(data, id.vars=names(data)[1:9], measure.vars=names(data)[10:59],
                                value.name="speed", variable.name="time")
  data_long[, time := as.numeric(gsub("vt","",time))]

  ## identify fish data that were removed or have incomplete entries
  rmv <- data[rowSums(is.na(data)) >0, names(data)[1:9], with=F]

  # fit poly2 to data
  fits <- lapply(concs, function(cnc) {

    # remove incomplete data
    to.fit <- data_long[!rmv, on=names(rmv)][conc==cnc]

    # fit to both phases: Light and Dark
    fit.L <- lm(speed ~ poly(time,2), data = to.fit[time%in%c(11:30)])
    fit.D <- lm(speed ~ poly(time,2), data = to.fit[time%in%c(31:50)])

    list("fit.L" = fit.L, "fit.D" = fit.D)
  })

  names(fits) <- concs
  fits
})
#### end #####

#### produce data.table of coefficient values and corresponding CIs ####

# coefficient data.table
coeff <- lapply(fits, lapply, lapply, `[[`, "coefficients")
coeff.dt <- lapply(coeff, lapply, unlist) %>%
            lapply(., lapply, unlist) %>%
            lapply(., do.call, what = 'rbind') %>%
            do.call('rbind', .) %>%
            as.data.table(., keep.rownames = "conc")

## format coefficient data.table
data.table::setnames(coeff.dt, names(coeff.dt)[-1], c("B_0.L","B_1.L","B_2.L","B_0.D","B_1.D","B_2.D"))
coeff.dt[, `:=` (cpid=c(rep(chms[1],10),rep(chms[2],8),rep(chms[3],7)), conc=as.numeric(conc))]
coeff.dt1 <- melt(coeff.dt, id.vars = c("cpid", "conc"), variable.name = "coeff")
coeff.dt1[grepl("L", coeff), phase := "Light"]
coeff.dt1[grepl("D", coeff), phase := "Dark"]
coeff.dt1[, coeff := gsub("\\..*","",coeff)]

# find confidence intervals, format and bind to data
CIs <- lapply(fits, lapply, lapply, confint)
CIs.dt <- lapply(CIs, lapply, do.call, what='rbind') %>%
            lapply(., do.call, what="rbind") %>%
            do.call("rbind", .) %>%
            as.data.table(., keep.rownames = "coeff")
CIs.dt[coeff=="(Intercept)", coeff := "B_0"]
CIs.dt[coeff=="poly(time, 2)1", coeff := "B_1"]
CIs.dt[coeff=="poly(time, 2)2", coeff := "B_2"]
CIs.dt[, phase:= rep(c("Light","Dark"), each=3, times=25)]

concs <- lapply(fits, names)
CIs.dt[, conc := as.numeric(c(rep(concs[[1]],each=6), rep(concs[[2]],each=6), rep(concs[[3]],each=6)))]
CIs.dt[, cpid := c(rep(rep(chms[1],times=length(concs[[1]])), each=6),
                   rep(rep(chms[2],times=length(concs[[2]])), each=6),
                   rep(rep(chms[3],times=length(concs[[3]])), each=6))]
data.table::setnames(CIs.dt, c("2.5 %","97.5 %"), c("min","max"))

## bind to coefficient data
coeff.dt2 <- coeff.dt1[CIs.dt, on = c("cpid","conc","phase","coeff")]
#### end ####
#### produce row vectors for tcplfitting ####

# start with Loperamide
chm <- chms[1]
lop_rows <- lapply(c("Light","Dark"), function(phs) {
              lapply(c("B_1","B_2"), function(param){
                data <- coeff.dt2[cpid==chm & coeff==param & phase==phs]
                bmed <- data[conc==0, value]
                resp <- data[conc!=0, value] - bmed
                conc <- data[conc!=0, conc]
                acid <- paste(param,phs,sep=".")
                cpid <- chm

                list(cpid = cpid,
                     acid = acid,
                     resp = resp,
                     conc = conc,
                     bmed = bmed,
                     data = data)
                })
              })
# 5-Fluorouracil
chm <- chms[2]
urcl_rows <- lapply(c("Light","Dark"), function(phs) {
  lapply(c("B_1","B_2"), function(param){
    data <- coeff.dt2[cpid==chm & coeff==param & phase==phs]
    bmed <- data[conc==0, value]
    resp <- data[conc!=0, value] - bmed
    conc <- data[conc!=0, conc]
    acid <- paste(param,phs,sep=".")
    cpid <- chm

    list(cpid = cpid,
         acid = acid,
         resp = resp,
         conc = conc,
         bmed = bmed,
         data = data)
  })
})
# 5,5-diphenylhydandoin
chm <- chms[3]
diphn_rows <- lapply(c("Light","Dark"), function(phs) {
  lapply(c("B_1","B_2"), function(param){
    data <- coeff.dt2[cpid==chm & coeff==param & phase==phs]
    bmed <- data[conc==0, value]
    resp <- data[conc!=0, value] - bmed
    conc <- data[conc!=0, conc]
    acid <- paste(param,phs,sep=".")
    cpid <- chm

    list(cpid = cpid,
         acid = acid,
         resp = resp,
         conc = conc,
         bmed = bmed,
         data = data)
  })
})

# combine
rows <- list("Loperamide" = lop_rows,
             "5-Fluorouracil" = urcl_rows,
             "5,5-diphenylhydandoin" = diphn_rows)
#### end ####
