#------------------------------------
# sudo tcplfit2 functons for fitting
#   Zachary Rowson
#   Rowson.Zachary@epa.gov
#   Created: 01/31/2022
#   Last Edit: 01/31/2022
#------------------------------------
#############################################################################################################
crc.local <- function(row,
                           fitmodels = c(
                             "cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
                             "exp4", "exp5"
                           ),
                           conthits = TRUE,
                           aicc = FALSE,
                           force.fit = TRUE, #EDITED
                           bidirectional = TRUE,
                           verbose = FALSE,
                           verbose.plot = TRUE,
                           do.plot = FALSE,
                           return.details = FALSE,
                           bmr_scale = 1.349,
                           bmd_low_bnd = NULL,
                           bmd_up_bnd = NULL) {
  # variable binding to pass cmd checks
  bmed <- cutoff <- onesd <- plot <- NULL
  # row needs to include cutoff and bmed
  # unpack row into the local environment, for ease: sample_id, dtxsid, casrn, name, time, pathway, size, conc, resp
  list2env(row, envir = environment())
  resp <- unlist(resp)
  conc <- unlist(conc)

  # prepare input EDIT: removed centralization about bmed, this is done is as_row
  conc <- conc[!is.na(resp)]
  resp <- resp[!is.na(resp)]
  identifiers <- row[!names(row) %in% c("conc", "resp", "bmed")] #EDIT

  # EDIT: calculate response medians in concRespCore to avoid running it in tcplfit2_coreZR and tcplggplotter
  logc <- log10(conc)
  rmds <- tapply(resp, logc, mean)

  # run the fits
  params <- tcplfit2_core.local(conc, resp, rmds,
                                  force.fit = conthits, bidirectional = bidirectional,
                                  fitmodels = fitmodels, verbose = verbose,
                                  do.plot = do.plot
  )

  # calculate the hitcall
  summary <- tcplhit2_core.local(params, conc, resp, bmr_scale, bmed, conthits, aicc, identifiers, bmd_low_bnd, bmd_up_bnd)

  # EDIT: create plotting summary if can.plot == TRUE
    plot <- tcplggplotter.local(resp, conc, row, rmds, bmed, params, summary, verbose.plot)

  if (return.details) {
    return(list(summary = summary, all.models = params, plot = plot))
  } else {
    return(list(summary = summary, plot = plot))
  }
}

#############################################################################################################################
tcplfit2_core.local <- function(conc, resp, rmds, force.fit = FALSE, bidirectional = TRUE, verbose = FALSE, do.plot = FALSE,
                            fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3", "exp4", "exp5"),
                            ...) {
  fitmodels <- unique(c("cnst", fitmodels)) # cnst models must be present for conthits but not chosen

  # first decide which of possible models will be fit
  modelnames <- c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3", "exp4", "exp5")
  # check for edge case where all responses are equal
  if(max(resp)==min(resp) && resp[1]==0){ # check if all response values are zero
    warning(paste("all response values are 0: add epsilon (1e-6) to all response elements.",
                  paste("\tResponse Range:",paste(range(resp),collapse = ",")),
                  sep = "\n")) # return a warning
    resp <- resp+1e-6 # adding epsilon to resp vector
  }
  # decide whether to run each model, then use generic functions to run model by name
  for (model in modelnames) {
    # only fit when four or more concentrations, the model is in fitmodels, and
    # ( either one response is above cutoff OR force.fit == T OR it's the constant model.)
    to.fit <- (length(rmds) >= 4 && model %in% fitmodels || force.fit || model == "cnst")
    fname <- paste0("fit", model) # requires each model function have name "fit____" where ____ is the model name
    # use do.call to call fit function; cnst has different inputs than others.
    assign(model, do.call(fname, list(
      conc = conc, resp = resp, bidirectional = bidirectional, verbose = verbose,
      nofit = !to.fit
    )))
    if (to.fit) {
      if (model %in% c("poly1", "poly2", "pow", "exp2", "exp3")) {
        # methods that grow without bound: top defined as model value at max conc
        assign(model, append(get(model), list(top = get(model)$modl[which.max(abs(get(model)$modl))]))) # top is taken to be highest model value
        assign(model, append(get(model), list(ac50 = acy(.5 * get(model)$top, get(model), type = model))))
      } else if (model %in% c("hill", "exp4", "exp5")) {
        # methods with a theoretical top/ac50
        assign(model, append(get(model), list(top = get(model)$tp)))
        assign(model, append(get(model), list(ac50 = get(model)$ga)))
      } else if (model == "gnls") {
        # gnls methods; use calculated top/ac50, etc.
        assign(model, append(get(model), list(top = acy(0, get(model), type = model, returntop = T))))
        # check if the theoretical top was calculated
        if(is.na(get(model)$top)){
          # if the theoretical top is NA return NA for ac50 and ac50_loss
          if(verbose){
            warning("'top' for 'gnls' is not able to be calculated returning NA.  AC50 for gain and loss directions are returned as NA.")
          }
          assign(model,append(get(model), list(ac50 = NA_real_,ac50_loss = NA_real_)))
        } else {
          assign(model, append(get(model), list(ac50 = acy(.5 * get(model)$top, get(model), type = model))))
          assign(model, append(get(model), list(ac50_loss = acy(.5 * get(model)$top, get(model), type = model, getloss = T))))
        }
      }
    }
    assign(model, append(get(model), list(func = gabi::tcplfit2_funcfitZR(model, fit=get(model))))) # EDIT attach fitted function with parameters to model for future curve plotting
  }
  # optionally print out AICs
  if (verbose) {
    print("aic values:")
    aics <- sapply(modelnames, function(x) {
      get(x)[["aic"]]
    })
    names(aics) <- modelnames
    print(aics)
    cat("Winner: ", modelnames[which.min(aics)])
  }

  # Produce logical object can.plot for tcplggplotter()
  shortnames <- modelnames[modelnames != "cnst"]
  successes <- sapply(shortnames, function(x) {
    get(x)[["success"]]
  })
  # EDIT below: declare can.plot object and append to out object
  if (do.plot && sum(successes, na.rm = T) == length(shortnames)) {
    can.plot <- TRUE
  } else {can.plot <- FALSE}

  # put all the model outputs into one list and return
  out <- c(
    mget(modelnames),
    list(can.plot = can.plot, modelnames = modelnames, ...)
  )

  return(out)
}

#############################################################################################################################

tcplhit2_core.local <- function(params, conc, resp, cutoff, onesd,bmr_scale = 1.349, bmed = 0, conthits = TRUE, aicc = FALSE,
                                identifiers = NULL, bmd_low_bnd = NULL, bmd_up_bnd = NULL) {
  # initialize parameters to NA
  a <- b <- tp <- p <- q <- ga <- la <- er <- top <- ac50 <- ac50_loss <- ac5 <- ac10 <- ac20 <- acc <- ac1sd <- bmd <- NA_real_
  bmdl <- bmdu <- caikwt <- mll <- NA_real_
  # get aics and degrees of freedom
  aics <- sapply(params$modelnames, function(x) {
    params[[x]][["aic"]]
  })
  dfs <- sapply(params$modelnames, function(x) {
    length(params[[x]][["pars"]])
  })
  aics <- aics[!is.na(aics)]
  if (sum(!is.na(aics)) == 0) {
    # if all fits failed, use none for method
    fit_method <- "none"
    rmse <- NA_real_
  } else {
    # use nested chisq to choose between poly1 and poly2, remove poly2 if it fails.
    # pvalue hardcoded to .05
    aics <- nestselect(aics, "poly1", "poly2", dfdiff = 1, pval = .05)
    dfs <- dfs[names(dfs) %in% names(aics)]
    # it's useful to keep original aics so we can extract loglikelihoods for nested models (above) and mll (below)
    aicc <- FALSE
    if (aicc) saics <- aics + 2 * dfs * (dfs + 1) / (length(resp) - dfs - 1) else saics <- aics
    conthits <- TRUE
    if (conthits) {
      # if all fits, except the constant fail, use none for the fit method
      # when continuous hit calling is in use
      if(sum(!is.na(aics)) == 1 & "cnst" %in% names(aics[!is.na(aics)])){
        fit_method <- "none"
        rmse <- NA_real_
      }else{
        # get aikaike weight of winner (vs constant) for cont hitcalls
        # never choose constant as winner for cont hitcalls
        nocnstaics <- saics[names(saics) != "cnst"]
        fit_method <- names(nocnstaics)[which.min(nocnstaics)]
        caikwt <- exp(-saics["cnst"] / 2) / (exp(-saics["cnst"] / 2) + exp(-saics[fit_method] / 2))
        if (is.nan(caikwt)) {
          term <- exp(saics["cnst"] / 2 - saics[fit_method] / 2)
          if (term == Inf) {
            caikwt <- 0
          } else {
            caikwt <- 1 / (1 + term)
          }
          # caikwt <- 1
        }
      }
    } else {
      fit_method <- names(saics)[which.min(saics)]
    }
    # if the fit_method is not reported as 'none' the obtain model information
    if(fit_method!="none"){
      fitout <- params[[fit_method]]
      rmse <- fitout$rme
      modpars <- fitout[fitout$pars]
      list2env(fitout, envir = environment()) # put all parameters in environment
    }
  }



  name.list <- c(
   "fit_method","rmse", "a", "b", "tp", "p", "q", "ga", "la", "er", "caikwt",
    "mll", "top")
  row <- as.data.frame(c(identifiers, mget(name.list)), stringsAsFactors = FALSE)
  return(row)
}

############################################################################################################################
tcplggplotter.local <- function(resp, conc, row, rmds, bmed, params, summary,
                                verbose.plot, unit.conc = paste0("\U03BC","M")) {

  data <- data.table::copy(row[["data"]])
  # plot concentration response curves and other statistics related to model of chemical activity

  # consolidate x & y coordinates for resp and bresp

  ## create data.table for responses and a pseudo concentration as x coordinate for bresp
  logc <- log10(conc[order(conc)])
  logc.sudo <- floor(min(logc)) - 1
  conc.sudo <- 10^(logc.sudo)
  concr <- c(conc[order(conc)], conc.sudo)
  logc <- c(logc, logc.sudo)

  # create x and y coordinates for plotting curve fits

  ## import curve fits
  list2env(params, env = environment())

  ## use winning fit's function to find f(x)
  fit_method <- summary$fit_method
  fit.func <- get(fit_method)[["func"]]
  x <- seq(from = conc.sudo, to = max(conc), by = 0.01)
  fit.xy <- data.table(fit = fit_method, conc = x, resp.hat = fit.func(x))

  # gather descriptive and inferential statistics for plotting

  ## round summary statistics for captions
  stats <- c("rmse", "caikwt", "top")
  cap <- round(as.numeric(summary[stats]), 3)
  names(cap) <- stats

  ## create caption
  if (verbose.plot) {
    caption <- paste0("Winning Fit=", fit_method,
                      ", top=", cap["top"],
                      ", rmse=", cap["rmse"],
                      ", caikwt= ", cap["caikwt"]

    )
  } else caption <- NULL

  ## create x-axis breaks and labels
  x.max <- floor(log10(max(concr)))
  x.breaks <- 10^(seq(from=min(logc), x.max))
  x.labels <- c("Control", format(x.breaks[-1], scientific=TRUE))

  ## create plot
  data[, `:=` (value = value-bmed,
                        min = min-bmed,
                        max = max-bmed)]
  data[conc==0, conc := conc.sudo]
  plot <- ggplot() +
    geom_point(data, mapping = aes(x=conc, y=value)) +
    geom_errorbar(data, mapping = aes(x=conc, ymin=min, ymax=max), width = 0.095) +
    geom_line(fit.xy, mapping = aes(x=conc, y=resp.hat, linetype="Winning Fit")) +
    labs(title = paste0(row$cpid, " for ", row$acid),
         x = unit.conc, # bquote("log"[10]~.(unit.conc)),
         y = expression(paste("Response - Vehicle Control Response")),
         caption = caption) +
    scale_x_continuous(trans = "log10", breaks = x.breaks, labels = x.labels) +
    theme_bw()

  return(plot)
}

# c("black", "cyan", "dark magenta", "red", "darkgoldenrod1",
#   "hotpink", "chartreuse", "darkred", "blue1")) +
