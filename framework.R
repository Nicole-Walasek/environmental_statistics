#this script contains the framework functions

# creates the framework

# helper function to compute to change points
# penalties can be AIC or MBIC, where MBIC applies a strciter change point criterion
# resulting in fewer changepoints

element_textbox_highlight <- function(..., hi.labels = NULL, hi.fill = NULL,
                                      hi.col = NULL, hi.box.col = NULL, hi.family = NULL) {
  structure(
    c(element_textbox(...),
      list(hi.labels = hi.labels, hi.fill = hi.fill, hi.col = hi.col, hi.box.col = hi.box.col, hi.family = hi.family)
    ),
    class = c("element_textbox_highlight", "element_textbox", "element_text", "element")
  )
}

element_grob.element_textbox_highlight <- function(element, label = "", ...) {
  if (label %in% element$hi.labels) {
    element$fill <- element$hi.fill %||% element$fill
    element$colour <- element$hi.col %||% element$colour
    element$box.colour <- element$hi.box.col %||% element$box.colour
    element$family <- element$hi.family %||% element$family
  }
  NextMethod()
}


# difference a sequence of numbers
# helper function
myDiff2 <- function(y, degree) {
  if (degree > 0) {
    for (idx in 1:degree) {
      y <- diff(y)
    }
  }
  return(y)
}

# extract random noise from TS
extractNoise <- function(y, freq) {
  dataTS <- decompose(ts(y, frequency = freq))
  return(dataTS$random)
  
}



# computes augmented dickey-fuller p values and outputs stationary TS
makeStationary <- function(y,
                           diffDegree = 1,
                           removeSeason = TRUE,
                           logTransform = FALSE,
                           type,
                           freq,
                           returnType) {
  
  data_SeasonalaAdj <- y
  if (logTransform) {
    data_SeasonalaAdj <- log(data_SeasonalaAdj)
  }
  
  
  dataTS <- decompose(ts(data_SeasonalaAdj, frequency = freq))
  if (removeSeason) {
    if (type == "additive") {
      data_SeasonalaAdj <-
        as.numeric(dataTS$x) - as.numeric(dataTS$seasonal)
      
    } else{
      data_SeasonalaAdj <-
        as.numeric(dataTS$x) / as.numeric(dataTS$seasonal)
    }
  } 
  
  data_diffAdj <- myDiff2(data_SeasonalaAdj, degree = diffDegree)
  
  if (returnType == 'dickey-fuller') {
    return(as.numeric(adf.test(data_diffAdj)$p.value))
  } else{
    return(data_diffAdj)
  }
  
}



# compute intensity / chronicity score
calcChronicityScore <- function(mean, var, theoreticalMax, type) {
  countMean <- plyr::count(mean)
  countMean$result <-
    (countMean$freq / sum(countMean$freq)) * countMean$x
  maxIDX <- which(countMean$freq == max(countMean$freq))
  maxVar <- var[which(mean == countMean$x[maxIDX])[1]]
  
  if (type == "chronicity") {
    score <- (countMean$result[maxIDX]) / (sqrt(maxVar) * theoreticalMax)
  } else if (type == 'intensity') {
    score <- (countMean$result[maxIDX]) / theoreticalMax
  } else{
    score <- mean(mean / sqrt(var))
  }
  
  return(max(score))
  
}


#get mean and variance per segment
changePoints_SegVals <-
  function(data, changePointsObj, type, segType) {
    x <- 1:length(data$y)
    
    cpsMean <- changePointsObj$cpMEAN@cpts
    cpsMean <- cpsMean+ 1
    
    cpsVar <- changePointsObj$cpVAR@cpts
    cpsVar <- cpsVar + 1
    
    cpsAll <- changePointsObj$cpNP_ALL@cpts
    cpsAll <- cpsAll + 1
    
    
    
    
    lMean <-
      split(x, rep(seq_along(cpsMean), times = diff(c(1, cpsMean))))
    lVar <- split(x, rep(seq_along(cpsVar), times = diff(c(1, cpsVar))))
    
    lAll <- split(x, rep(seq_along(cpsAll), times = diff(c(1, cpsAll))))
    
    
    if (type == "All") {
      resMean <- sapply(lAll, function(x)
        mean(data$y[x]))
      resVar <- sapply(lAll, function(x)
        sd(data$y[x]))
      resVarMean <- resMean
      cpsMean <- cpsAll
      cpsVar <- cpsAll
      
    } else{
      resMean <- sapply(lMean, function(x)
        mean(data$y[x]))
      resVar <- sapply(lVar, function(x)
        sd(data$y[x]))
      resVarMean <- sapply(lVar, function(x)
        mean(data$y[x]))
    }
    
    resultMean <- as.numeric(rep(resMean, times = diff(c(1, cpsMean))))
    resultVar <- as.numeric(rep(resVar, times = diff(c(1, cpsVar))))
    resultVarMean <-
      as.numeric(rep(resVarMean, times = diff(c(1, cpsVar))))
    if (segType == 'Mean') {
      return(resultMean)
    } else if (segType == 'Var') {
      return(resultVar)
    } else{
      return(resultVarMean)
    }
    
  }


#get mean and variance per segment
changePoints_SegValsVar <-
  function(data, changePointsObj, type, segType) {
    x <- 1:length(data$yVar)
    
    cpsMean <- changePointsObj$cpMEAN@cpts
    cpsMean <- cpsMean+ 1
    
    cpsVar <- changePointsObj$cpVAR@cpts
    cpsVar <- cpsVar + 1
    
    cpsAll <- changePointsObj$cpNP_ALL@cpts
    cpsAll <- cpsAll + 1
    
    
    
    
    lMean <-
      split(x, rep(seq_along(cpsMean), times = diff(c(1, cpsMean))))
    lVar <- split(x, rep(seq_along(cpsVar), times = diff(c(1, cpsVar))))
    
    lAll <- split(x, rep(seq_along(cpsAll), times = diff(c(1, cpsAll))))
    
    
    if (type == "All") {
      resMean <- sapply(lAll, function(x)
        mean(data$yVar[x]))
      resVar <- sapply(lAll, function(x)
        sd(data$yVar[x]))
      resVarMean <- resMean
      cpsMean <- cpsAll
      cpsVar <- cpsAll
      
    } else{
      resMean <- sapply(lMean, function(x)
        mean(data$yVar[x]))
      resVar <- sapply(lVar, function(x)
        sd(data$yVar[x]))
      resVarMean <- sapply(lVar, function(x)
        mean(data$yVar[x]))
    }
    
    resultMean <- as.numeric(rep(resMean, times = diff(c(1, cpsMean))))
    resultVar <- as.numeric(rep(resVar, times = diff(c(1, cpsVar))))
    resultVarMean <-
      as.numeric(rep(resVarMean, times = diff(c(1, cpsVar))))
    if (segType == 'Mean') {
      return(resultMean)
    } else if (segType == 'Var') {
      return(resultVar)
    } else{
      return(resultVarMean)
    }
    
  }


identifyChangePoints <- function(data, penalty, minseglen) {
  cpMEAN <-
    cpt.mean(data$y,
             method = 'PELT',
             penalty = penalty,
             minseglen = minseglen)
  cpVAR <-
    cpt.var(data$y,
            method = 'PELT',
            penalty = penalty,
            minseglen = minseglen)
  cpNP_ALL <-
    cpt.np(data$y,
           method = 'PELT',
           penalty = penalty,
           minseglen = minseglen)
  result <-
    list(
      "cpMEAN" = cpMEAN,
      "cpVAR" = cpVAR,
      "cpNP_ALL" = cpNP_ALL,
      minseglen = minseglen
    )
  return(result)
  
}

pos_neg_changePoints <-
  function(data, changePointsObj, type, direction) {
    x <- 1:length(data$y)
    if (type == 'mean') {
      cps <- changePointsObj$cpMEAN@cpts
    } else{
      cps <- changePointsObj$cpVAR@cpts
    }
    cps[length(cps)] <- cps[length(cps)] + 1
    
    
    pat <- pat <- rep(seq_along(cps), times = diff(c(1, cps)))
    l <- split(x, pat)
    
    if (type == 'mean') {
      res <- sapply(l, function(x)
        mean(data$y[x]))
      
    } else{
      res <- sapply(l, function(x)
        var(data$y[x]))
      
    }
    diffRes <- diff(res)
    
    if (direction == 'positive') {
      return(sum(diffRes > 0))
    } else{
      return(sum(diffRes < 0))
    }
    
  }

# also add the average distance between change points and the corresponding standard deviation

distance_changePoints <-
  function(data, changePointsObj, type, stats) {
    x <- 1:length(data$y)
    if (type == 'mean') {
      cps <- changePointsObj$cpMEAN@cpts
    } else if (type == 'var') {
      cps <- changePointsObj$cpVAR@cpts
    } else{
      cps <- changePointsObj$cpNP_ALL@cpts
    }
    
    if (length(cps) == 1) {
      if (stats == 'mean' | stats == 'max') {
        return(cps[1])
      } else{
        return(0)
      }
    }
    cps <- c(1, cps)
    diffCps <- diff(cps)
    if (stats == 'mean') {
      return(mean(diffCps))
    } else if (stats == 'sd') {
      return(sd(diffCps))
    } else{
      return(max(diffCps))
    }
  }

spectral_analysis <- function(data, freq) {
  data.spect <-
    spectrum(na.omit(extractNoise(data$y, freq = freq)), log = "no", plot = FALSE)
  data.spect$freq_log <- log(data.spect$freq)
  data.spect$spec_log <- log(data.spect$spec)
  return(data.spect)
}

spectralIdentifyCycle <- function(dataRaw) {
  # return be extended to return the n biggest cycles
  data <- spectrum(dataRaw$y, log = "no", plot = FALSE)
  cycle <- 1 / data$freq[which.max(data$spec)]
  if (cycle == data$n.used) {
    return(0)
  } else{
    return(cycle)
  }
  
}

spectralModel <- function(data) {
  #specify the model of the mean here
  lm(data$spec_log ~ data$freq_log)
}

extraxSpectralCoef <- function(data) {
  return(data$coefficients[2]*-1)
}


framework.applyModel <- function(data,
                                 acLAG,
                                 minseglen,
                                 percentile,
                                 meanModel = function(data) {
                                   lm(y ~ time, data = data)
                                 },
                                 varModel = function(data) {
                                   lm(yVar ~ time, data = data)
                                 },
                                 residuals = FALSE,
                                 sqaured = FALSE,
                                 diffDegree = 1,
                                 removeSeason = TRUE,
                                 logTransform = FALSE,
                                 freq,
                                 type) {
  data_nested <- data %>%
    group_by(ppID) %>%
    nest()
  
  data_models <- data_nested %>%
    group_by(ppID) %>%
    mutate(
      #mean stats
      data_lm      = map(data, meanModel),
      data_stationary = map(data, function(data) {
        makeStationary(
          data$y,
          diffDegree = diffDegree,
          removeSeason = removeSeason,
          logTransform = logTransform,
          type = type,
          freq = freq,
          returnType = 'ts'
        )
      }),
      
      data_adf = map_dbl(data, function(data) {
        makeStationary(
          data$y,
          diffDegree = 0,
          removeSeason = FALSE,
          logTransform = FALSE,
          type = type,
          freq = freq,
          returnType = 'dickey-fuller'
        )
      }),
      
      data_stationary_adf = map_dbl(data, function(data) {
        makeStationary(
          data$y,
          diffDegree = diffDegree,
          removeSeason = removeSeason,
          logTransform = logTransform,
          type = type,
          freq = freq,
          returnType = 'dickey-fuller'
        )
      }),
      
      data_acStat      = map(data_stationary, function(data) {
        result <-
          data.frame(
            acf(
              data,
              lag.max = length(data) - 1,
              na.action = na.pass,
              plot = FALSE
            )$acf[2:(length(data))],
            pacf(
              data,
              lag.max = length(data) - 1,
              na.action = na.pass,
              plot = FALSE
            )$acf[1:(length(data) - 1)]
          )
        
        names(result) <- c("acf", "pacf")
        result$lag <- as.factor(c(1:(length(data) - 1)))
        result <- as_tibble(result)
        
      }),
      data_ac      = map(data, function(data) {
        result <-
          data.frame(
            acf(
              as.numeric(data$y),
              lag.max = nrow(data) - 1,
              plot = FALSE
            )$acf[2:(nrow(data))],
            pacf(as.numeric(data$y), plot = FALSE)$acf[1:(nrow(data) -
                                                            1)]
          )
        names(result) <- c("acf", "pacf")
        result$lag <- as.factor(c(1:(nrow(data) - 1)))
        result <- as_tibble(result)
      }),
      lm_statsMean = map(data_lm, glance),
      lm_coefsMean = map(data_lm, tidy),
      lm_valsMean  = map(data_lm, augment)
    )
  
  # add variance in the outcome variable
  if (residuals) {
    data <- ddply(data, .(ppID)) %>%
      mutate(yVar = unlist(map(data_models$lm_valsMean, ".resid")))
    
  } else{
    data <- ddply(data, .(ppID)) %>%
      mutate(yVar =
               as.numeric(sapply(split(
                 data$y, data$ppID
               ), function(x) {
                 as.numeric(scale(x, scale = FALSE, center = TRUE))
               })))
  }
  
  if (sqaured) {
    data$yVar <- data$yVar ** 2
  }
  
  # adding the residuals
  data <- data %>%
    mutate(yRes = unlist(map(data_models$lm_valsMean, ".resid")))
  
  data_nested <- data %>%
    group_by(ppID) %>%
    nest()
  data_models$data <- data_nested$data
  
  
  data_models <- data_models %>%
    group_by(ppID) %>%
    mutate(
      # variance statistics
      dataVar_lm      = map(data, varModel),
      
      
      dataVar_stationary = map(data, function(data) {
        makeStationary(
          data$yVar,
          diffDegree = diffDegree,
          removeSeason = removeSeason,
          logTransform = logTransform,
          type = type,
          freq = freq,
          returnType = 'ts'
        )
      }),
      
      dataVar_adf = map_dbl(data, function(data) {
        makeStationary(
          data$yVar,
          diffDegree = 0,
          removeSeason = FALSE,
          logTransform = FALSE,
          type = type,
          freq = freq,
          returnType = 'dickey-fuller'
        )
      }),
      
      dataVar_stationary_adf = map_dbl(data, function(data) {
        makeStationary(
          data$yVar,
          diffDegree = diffDegree,
          removeSeason = removeSeason,
          logTransform = logTransform,
          type = type,
          freq = freq,
          returnType = 'dickey-fuller'
        )
      }),
      
      dataVar_acStat      = map(dataVar_stationary, function(data) {
        result <-
          data.frame(
            acf(
              data,
              lag.max = length(data) - 1,
              na.action = na.pass,
              plot = FALSE
            )$acf[2:(length(data))],
            pacf(
              data,
              lag.max = length(data) - 1,
              na.action = na.pass,
              plot = FALSE
            )$acf[1:(length(data) - 1)]
          )
        
        names(result) <- c("acf", "pacf")
        result$lag <- as.factor(c(1:(length(data) - 1)))
        result <- as_tibble(result)
      }),
      
      dataVar_ac      = map(data, function(data) {
        result <-
          data.frame(
            acf(
              as.numeric(data$yVar),
              lag.max = nrow(data) - 1,
              plot = FALSE
            )$acf[2:(nrow(data))],
            pacf(as.numeric(data$yVar), plot = FALSE)$acf[1:(nrow(data) -
                                                               1)]
          )
        names(result) <- c("acf", "pacf")
        result$lag <- as.factor(c(1:(nrow(data) - 1)))
        result <- as_tibble(result)
      }),
      
      lm_statsVar = map(dataVar_lm, glance),
      lm_coefsVar = map(dataVar_lm, tidy),
      lm_valsVar  = map(dataVar_lm, augment),
      
      # add changepoints
      changePoints_MBIC = pmap(list(data, "MBIC", minseglen), identifyChangePoints),
      changePoints_AIC  = pmap(list(data, "AIC", minseglen), identifyChangePoints),
      
      # add spectral analysis of the noise in the data
      dataSpectral = pmap(list(data, freq), spectral_analysis),
      dataSpectralLM = map(dataSpectral, spectralModel),
      dataMaxPeriod = map_dbl(data, spectralIdentifyCycle),
      
      
    )
  
  #adding all the other model coefficients for the mean
  
  # add model coefficients
  
  #add ppID
  ppID <- data_models$ppID
  df <- data.frame(ppID)
  df$m_intercept <-
    map_dbl(data_models$lm_coefsMean, function(x)
      unlist(x[1, "estimate"]))
  degree = nrow(data_models$lm_coefsMean[[1]]) - 1
  for (idx in 1:degree) {
    current <-
      map_dbl(data_models$lm_coefsMean, function(x)
        unlist(x[idx + 1, "estimate"]))
    df <- cbind(df, current)
    
    names(df)[names(df) == "current"] <-  sprintf("m_slope%s", idx)
  }
  
  
  #add rsquared
  df$m_rsqrd <- data_models$lm_statsMean %>% map_dbl("r.squared")
  
  # add mean and sd
  df$m_mean <-
    map_dbl(data_models$data, function(data) {
      mean(data$y)
    })
  df$m_sd <- map_dbl(data_models$data, function(data) {
    sd(data$y)
  })
  
  #add max, min, range, IQR
  df$m_max <- map_dbl(data_models$data, function(data) {
    max(data$y)
  })
  df$m_min <- map_dbl(data_models$data, function(data) {
    min(data$y)
  })
  df$m_range <-
    map_dbl(data_models$data, function(data) {
      max(data$y) - min(data$y)
    })
  df$m_IQR <-
    map_dbl(data_models$data, function(data) {
      IQR(data$y)
    })
  
  # add approximate entropy
  df$m_apprEntropy <-
    map_dbl(data_models$data, function(data) {
      approx_entropy(data$y)
    })
  
  # add spectral coefficient THIS NEEDS TO CHANGE
  df$m_spectralCoef = map_dbl(data_models$dataSpectralLM, extraxSpectralCoef)
  
  
  dfTibble <- as_tibble(df)
  # adding all the ac and pac values
  dfAC <-
    data.frame(data.frame(matrix(
      ncol = 0, nrow = nrow(data_models)
    )))
  dfPAC <-
    data.frame(data.frame(matrix(
      ncol = 0, nrow = nrow(data_models)
    )))
  for (idx in 1:acLAG) {
    current <-
      map_dbl(data_models$data_ac, function(x)
        unlist(x[idx, "acf"]))
    dfAC <- cbind(dfAC, current)
    names(dfAC)[names(dfAC) == "current"] <-  sprintf("m_ac%s", idx)
    
    current <-
      map_dbl(data_models$data_ac, function(x)
        unlist(x[idx, "pacf"]))
    dfPAC <- cbind(dfPAC, current)
    names(dfPAC)[names(dfPAC) == "current"] <-
      sprintf("m_pac%s", idx)
  }
  dfACTibble <- as_tibble(dfAC)
  dfPACTibble <- as_tibble(dfPAC)
  meanStats = bind_cols("model_coef" = dfTibble,
                        "AC" = dfACTibble,
                        "PAC" = dfPACTibble)
  
  
  # statistics of the stationary distribution
  #add ppID
  ppID <- data_models$ppID
  df <- data.frame(ppID)
  # add mean and sd
  df$stat_mean <-
    map_dbl(data_models$data_stationary, function(data) {
      mean(data)
    })
  df$stat_sd <-
    map_dbl(data_models$data_stationary, function(data) {
      sd(data)
    })
  
  #add max, min, range, IQR
  df$stat_max <-
    map_dbl(data_models$data_stationary, function(data) {
      max(data)
    })
  df$stat_min <-
    map_dbl(data_models$data_stationary, function(data) {
      min(data)
    })
  df$stat_range <-
    map_dbl(data_models$data_stationary, function(data) {
      max(data) - min(data)
    })
  df$stat_IQR <-
    map_dbl(data_models$data_stationary, function(data) {
      IQR(data)
    })
  # add approximate entropy
  df$stat_apprEntropy <-
    map_dbl(data_models$data_stationary, function(data) {
      approx_entropy(data)
    })
  
  # all of the above for stationary dist of the variance
  # add mean and sd
  df$statVar_mean <-
    map_dbl(data_models$dataVar_stationary, function(data) {
      mean(data)
    })
  df$statVar_sd <-
    map_dbl(data_models$dataVar_stationary, function(data) {
      sd(data)
    })
  
  #add max, min, range, IQR
  df$statVar_max <-
    map_dbl(data_models$dataVar_stationary, function(data) {
      max(data)
    })
  df$statVar_min <-
    map_dbl(data_models$dataVar_stationary, function(data) {
      min(data)
    })
  df$statVar_range <-
    map_dbl(data_models$dataVar_stationary, function(data) {
      max(data) - min(data)
    })
  df$statVar_IQR <-
    map_dbl(data_models$dataVar_stationary, function(data) {
      IQR(data)
    })
  # add approximate entropy
  df$stat_apprEntropy <-
    map_dbl(data_models$dataVar_stationary, function(data) {
      approx_entropy(data)
    })
  
  
  
  dfTibble <- as_tibble(df)
  
  #stationary dist of the raw data
  dfACStat <-
    data.frame(data.frame(matrix(
      ncol = 0, nrow = nrow(data_models)
    )))
  dfPACStat <-
    data.frame(data.frame(matrix(
      ncol = 0, nrow = nrow(data_models)
    )))
  for (idx in 1:acLAG) {
    current <-
      map_dbl(data_models$data_acStat, function(x)
        unlist(x[idx, "acf"]))
    dfACStat <- cbind(dfACStat, current)
    names(dfACStat)[names(dfACStat) == "current"] <-
      sprintf("stat_ac%s", idx)
    
    current <-
      map_dbl(data_models$data_acStat, function(x)
        unlist(x[idx, "pacf"]))
    dfPACStat <- cbind(dfPACStat, current)
    names(dfPACStat)[names(dfPACStat) == "current"] <-
      sprintf("stat_pac%s", idx)
  }
  dfACStatTibble <- as_tibble(dfACStat)
  dfPACStatTibble <- as_tibble(dfPACStat)
  
  # stationary dist of the variance
  dfACStatVar <-
    data.frame(data.frame(matrix(
      ncol = 0, nrow = nrow(data_models)
    )))
  dfPACStatVar <-
    data.frame(data.frame(matrix(
      ncol = 0, nrow = nrow(data_models)
    )))
  for (idx in 1:acLAG) {
    current <-
      map_dbl(data_models$dataVar_acStat, function(x)
        unlist(x[idx, "acf"]))
    dfACStatVar <- cbind(dfACStatVar, current)
    names(dfACStatVar)[names(dfACStatVar) == "current"] <-
      sprintf("statVar_ac%s", idx)
    
    current <-
      map_dbl(data_models$dataVar_acStat, function(x)
        unlist(x[idx, "pacf"]))
    dfPACStatVar <- cbind(dfPACStatVar, current)
    names(dfPACStatVar)[names(dfPACStatVar) == "current"] <-
      sprintf("statVar_pac%s", idx)
  }
  dfACStatVarTibble <- as_tibble(dfACStatVar)
  dfPACStatVarTibble <- as_tibble(dfPACStatVar)
  
  
  
  
  
  stationaryStats = bind_cols(
    "stat_stats" = dfTibble,
    "AC_stat" = dfACStatTibble,
    "PAC_Stat" = dfPACStatTibble,
    "AC_StatVar" = dfACStatVarTibble,
    "PAC_StatVar" = dfPACStatVarTibble
  )
  
  #adding all the other model coefficients for the variance
  ppID <- data_models$ppID
  df <- data.frame(ppID)
  df$v_intercept <-
    map_dbl(data_models$lm_coefsVar, function(x)
      unlist(x[1, "estimate"]))
  degree = nrow(data_models$lm_coefsVar[[1]]) - 1
  for (idx in 1:degree) {
    current <-
      map_dbl(data_models$lm_coefsVar, function(x)
        unlist(x[idx + 1, "estimate"]))
    df <- cbind(df, current)
    
    names(df)[names(df) == "current"] <-  sprintf("v_slope%s", idx)
  }
  
  #add rsquared
  df$v_rsqrd <- data_models$lm_statsVar %>% map_dbl("r.squared")
  
  # add mean and sd
  df$v_mean <-
    map_dbl(data_models$data, function(data) {
      mean(data$yVar)
    })
  df$v_sd <-
    map_dbl(data_models$data, function(data) {
      sd(data$yVar)
    })
  
  #add max, min, range, IQR
  df$v_max <-
    map_dbl(data_models$data, function(data) {
      max(data$yVar)
    })
  df$v_min <-
    map_dbl(data_models$data, function(data) {
      min(data$yVar)
    })
  df$v_range <-
    map_dbl(data_models$data, function(data) {
      max(data$yVar) - min(data$yVar)
    })
  df$v_IQR <-
    map_dbl(data_models$data, function(data) {
      IQR(data$yVar)
    })
  
  
  # add approximate entropy
  df$v_apprEntropy <-
    map_dbl(data_models$data, function(data) {
      approx_entropy(data$yVar)
    })
  
  #add ppID
  df$ppID <- data_models$ppID
  
  dfTibble <- as_tibble(df)
  # adding all the ac values
  dfAC <-
    data.frame(data.frame(matrix(
      ncol = 0, nrow = nrow(data_models)
    )))
  dfPAC <-
    data.frame(data.frame(matrix(
      ncol = 0, nrow = nrow(data_models)
    )))
  for (idx in 1:acLAG) {
    current <-
      map_dbl(data_models$dataVar_ac, function(x)
        unlist(x[idx, "acf"]))
    dfAC <- cbind(dfAC, current)
    names(dfAC)[names(dfAC) == "current"] <-  sprintf("v_ac%s", idx)
    
    current <-
      map_dbl(data_models$dataVar_ac, function(x)
        unlist(x[idx, "pacf"]))
    dfPAC <- cbind(dfPAC, current)
    names(dfPAC)[names(dfPAC) == "current"] <-
      sprintf("v_pac%s", idx)
  }
  dfACTibble <- as_tibble(dfAC)
  dfPACTibble <- as_tibble(dfPAC)
  
  varStats = bind_cols("model_coef" = dfTibble,
                       "AC" = dfACTibble,
                       "PAC" = dfPACTibble)
  
  
  #adding all the other model coefficients for the residuals
  
  ppID <- data_models$ppID
  df <- data.frame(ppID)
  
  # add mean and sd
  df$r_mean <-
    map_dbl(data_models$data, function(data) {
      mean(data$yRes)
    })
  df$r_sd <-
    map_dbl(data_models$data, function(data) {
      sd(data$yRes)
    })
  
  #add max, min, range, IQR
  df$r_max <-
    map_dbl(data_models$data, function(data) {
      max(data$yRes)
    })
  df$r_min <-
    map_dbl(data_models$data, function(data) {
      min(data$yRes)
    })
  df$r_range <-
    map_dbl(data_models$data, function(data) {
      max(data$yRes) - min(data$yRes)
    })
  df$r_IQR <-
    map_dbl(data_models$data, function(data) {
      IQR(data$yRes)
    })
  
  
  # add approximate entropy
  df$r_apprEntropy <-
    map_dbl(data_models$data, function(data) {
      approx_entropy(data$yRes)
    })
  
  
  dfTibble <- as_tibble(df)
  
  resStats = bind_cols("res_stats" = dfTibble)
  
  # create change point statistics for MBIC
  meanNumberMBIC <-
    map_dbl(data_models$changePoints_MBIC, function(x) {
      length(x$cpMEAN@cpts) - 1
    })
  dfMBIC <- data.frame(meanNumberMBIC)
  dfMBIC$posMeanNumberMBIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        'mean',
        'positive'
      ),
      pos_neg_changePoints
    )
  dfMBIC$negMeanNumberMBIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        'mean',
        'negative'
      ),
      pos_neg_changePoints
    )
  
  dfMBIC$varNumberMBIC <-
    map_dbl(data_models$changePoints_MBIC, function(x) {
      length(x$cpVAR@cpts) - 1
    })
  dfMBIC$posVarNumberMBIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        'var',
        'positive'
      ),
      pos_neg_changePoints
    )
  dfMBIC$negVarNumberMBIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        'var',
        'negative'
      ),
      pos_neg_changePoints
    )
  
  dfMBIC$npNumberMBIC <-
    map_dbl(data_models$changePoints_MBIC, function(x) {
      length(x$cpNP_ALL@cpts) - 1
    })
  
  dfMBIC$avg_Distance_Mean_MBIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        'mean',
        'mean'
      ),
      distance_changePoints
    )
  dfMBIC$sd_Distance_Mean_MBIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        'mean',
        'sd'
      ),
      distance_changePoints
    )
  dfMBIC$max_Distance_Mean_MBIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        'mean',
        'max'
      ),
      distance_changePoints
    )
  dfMBIC$avg_Distance_Var_MBIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        'var',
        'mean'
      ),
      distance_changePoints
    )
  dfMBIC$sd_Distance_Var_MBIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        'var',
        'sd'
      ),
      distance_changePoints
    )
  dfMBIC$max_Distance_Var_MBIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        'var',
        'max'
      ),
      distance_changePoints
    )
  dfMBIC$avg_Distance_np_MBIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        'np',
        'mean'
      ),
      distance_changePoints
    )
  dfMBIC$sd_Distance_np_MBIC <-
    pmap_dbl(
      list(data_models$data, data_models$changePoints_MBIC, 'np', 'sd'),
      distance_changePoints
    )
  dfMBIC$max_Distance_np_MBIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        'np',
        'max'
      ),
      distance_changePoints
    )
  
  dfMBIC$ppID <- ppID
  dfMBIC <- as_tibble(dfMBIC)
  
  
  
  meanNumberAIC <-
    map_dbl(data_models$changePoints_AIC, function(x) {
      length(x$cpMEAN@cpts) - 1
    })
  dfAIC <- data.frame(meanNumberAIC)
  
  dfAIC$posMeanNumberAIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_AIC,
        'mean',
        'positive'
      ),
      pos_neg_changePoints
    )
  dfAIC$negMeanNumberAIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_AIC,
        'mean',
        'negative'
      ),
      pos_neg_changePoints
    )
  
  dfAIC$varNumberAIC <-
    map_dbl(data_models$changePoints_AIC, function(x) {
      length(x$cpVAR@cpts) - 1
    })
  dfAIC$posVarNumberAIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_AIC,
        'var',
        'positive'
      ),
      pos_neg_changePoints
    )
  dfAIC$negVarNumberAIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_AIC,
        'var',
        'negative'
      ),
      pos_neg_changePoints
    )
  
  dfAIC$npNumberAIC <-
    map_dbl(data_models$changePoints_AIC, function(x) {
      length(x$cpNP_ALL@cpts) - 1
    })
  
  dfAIC$avg_Distance_Mean_AIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_AIC,
        'mean',
        'mean'
      ),
      distance_changePoints
    )
  dfAIC$sd_Distance_Mean_AIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_AIC,
        'mean',
        'sd'
      ),
      distance_changePoints
    )
  dfAIC$max_Distance_Mean_AIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_AIC,
        'mean',
        'max'
      ),
      distance_changePoints
    )
  dfAIC$avg_Distance_Var_AIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_AIC,
        'var',
        'mean'
      ),
      distance_changePoints
    )
  dfAIC$sd_Distance_Var_AIC <-
    pmap_dbl(
      list(data_models$data, data_models$changePoints_AIC, 'var', 'sd'),
      distance_changePoints
    )
  dfAIC$max_Distance_Var_AIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_AIC,
        'var',
        'max'
      ),
      distance_changePoints
    )
  dfAIC$avg_Distance_np_AIC <-
    pmap_dbl(
      list(
        data_models$data,
        data_models$changePoints_AIC,
        'np',
        'mean'
      ),
      distance_changePoints
    )
  dfAIC$sd_Distance_np_AIC <-
    pmap_dbl(
      list(data_models$data, data_models$changePoints_AIC, 'np', 'sd'),
      distance_changePoints
    )
  dfAIC$max_Distance_np_AIC <-
    pmap_dbl(
      list(data_models$data, data_models$changePoints_AIC, 'np', 'max'),
      distance_changePoints
    )
  
  dfAIC$ppID <- ppID
  dfAIC <- as_tibble(dfAIC)
  
  # add chronicity and intensity score measures
  data <- data_models %>%
    dplyr::select(data, ppID) %>%
    unnest(cols = c(data))
  
  theoreticalMax <-
    as.numeric(quantile(data$y, probs = c(percentile)))
  
  cpVar <-
    pmap(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        "All",
        "Var"
      ),
      changePoints_SegVals
    )
  cpMean <-
    pmap(
      list(
        data_models$data,
        data_models$changePoints_MBIC,
        "All",
        "Mean"
      ),
      changePoints_SegVals
    )
  chronicity_score <-
    unlist(pmap(
      list(cpMean, cpVar, theoreticalMax, 'chronicity'),
      calcChronicityScore
    ))
  intensity_score <-
    unlist(pmap(
      list(cpMean, cpVar, theoreticalMax, 'intensity'),
      calcChronicityScore
    ))
  scaled_mean <-
    unlist(pmap(
      list(cpMean, cpVar, theoreticalMax, 'scaled mean'),
      calcChronicityScore
    ))
  compressed_statistics <- data.frame(chronicity_score)
  compressed_statistics$intensity_score <- intensity_score
  compressed_statistics$scaled_mean <- scaled_mean
  #could add more compressed statistics here
  compressed_statistics <- as_tibble(compressed_statistics)
  
  data_models <- data_models %>%
    add_column("Statistics_Mean" = meanStats)
  
  data_models <- data_models %>%
    add_column("Statistics_Variance" = varStats)
  
  data_models <- data_models %>%
    add_column("Statistics_Stationary" = stationaryStats)
  
  
  data_models <- data_models %>%
    add_column("Statistics_Residuals" = resStats)
  
  data_models <- data_models %>%
    add_column("Statistics_ChangePointsMBIC" = dfMBIC)
  
  data_models <- data_models %>%
    add_column("Statistics_ChangePointsAIC" = dfAIC)
  
  data_models <- data_models %>%
    add_column("Compressed_Statistics" = compressed_statistics)
  
  return(data_models)
  
}


# plots the framework models

framework.plotModel <-
  function(dataAll,
           numPP,
           modelTitle,
           type,
           freq = 1,
           xlab = "time",
           ylab,
           segment = TRUE,
           sampleList = c()) {
    
    timeDF <- dataAll %>%
      dplyr::select(data, ppID) %>%
      unnest(cols = c(data)) %>%
      dplyr::select(time, ppID)
    
    timePlot <- timeDF$time
    
    if (type == "mean") {
      # plot mean
      data <- dataAll %>%
        dplyr::select(lm_valsMean, ppID) %>%
        unnest(cols = c(lm_valsMean)) %>%
        add_column(timePlot = timePlot)
      
      data <- as.data.frame(data)
      
    } else{
      # plot variance
      data <- dataAll %>%
        dplyr::select(lm_valsVar, ppID) %>%
        unnest(cols = c(lm_valsVar)) %>%
        add_column(timePlot = timePlot)
      
      data <- as.data.frame(data)
      
      
    }
    
    if (length(sampleList) == 0) {
      results <-
        filter(data, ppID %in% sample(unique(data$ppID), numPP)) #select random participants for plotting
    } else{
      results <- filter(data, ppID %in% sampleList)
    }
    
    
    results$ppID <- as.factor(results$ppID)
    names(results)[names(results) == '.fitted'] <- 'fittedTimeModel'
    names(results)[names(results) == '.resid'] <-
      'residualsTimeModel'
    
    if (type == "mean") {
      names(results)[names(results) == 'y'] <- 'originalData'
      yLab <- ylab
    } else{
      names(results)[names(results) == 'yVar'] <- 'originalData'
      yLab <- sprintf("%s variance", ylab)
    }
    
    
    
    pModel <- ggplot(data = results)
    
    if (segment) {
      pModel <-
        pModel + geom_segment(
          aes(
            x = timePlot,
            xend = timePlot,
            y = originalData,
            yend = originalData - residualsTimeModel
          ),
          color = "grey"
        )
    }
    
    pModel <- pModel +
      geom_line(aes(x = timePlot, y = originalData, group = ppID),
                color = "black",
                size = 0.5) +
      geom_point(aes(x =
                       
                       timePlot, y = originalData, group = ppID), size = (1 /
                                                                            (0.7 * numPP))) + geom_line(
                                                                              aes(x = timePlot, y = fittedTimeModel, group = ppID),
                                                                              size = 0.8,
                                                                              color = "firebrick"
                                                                            )
    pModel <-
      pModel + labs(y = yLab, x = xlab) + scale_x_continuous(breaks =
                                                               seq(0, nObs, freq), labels = c(0:(nObs / freq)))
    
    pModel <- pModel + labs(title = paste(modelTitle, " model"))
    pModel <-
      pModel +       theme(legend.position = "none") + theme_bw() + theme(
        plot.title = element_text(face = "bold", size = 15) ,
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey70', size = 1),
        legend.position = c(0.06, 0.75),
        text = element_text(size = 15)
      )
    
    pRes <-
      ggplot(data = results, aes(x = timePlot, y = residualsTimeModel, group = ppID)) + geom_line(color = "#535754", aes(group = ppID), size = 0.5) + geom_point(aes(group = ppID), size = (1 /
                                                                                                                                                                                              (0.7 *
                                                                                                                                                                                                 numPP)))
    pRes <-
      pRes + labs(y = "residual", x = xlab) + scale_x_continuous(breaks =
                                                                   seq(0, nObs, freq), labels = c(0:(nObs / freq)))
    pRes <- pRes + labs(title = "residual error series")
    pRes <-
      pRes +       theme(legend.position = "none") + theme_bw() + theme(
        plot.title = element_text(face = "bold", size = 15) ,
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey70', size = 1),
        legend.position = c(0.06, 0.75),
        text = element_text(size = 15)
      )
    
    return(list('pModel' = pModel, 'pRes' = pRes))
  }


framework.plotColorOfNoise <- function(y, binwidth) {
  df_spectralCoef <- data.frame(y)
  names(df_spectralCoef) <- c("beta")
  color <- rep("", nrow(df_spectralCoef))
  color[df_spectralCoef < 0] <- "blue"
  color[df_spectralCoef > -0.5 & df_spectralCoef < 0.5] <- "white"
  color[df_spectralCoef >= 0.5 & df_spectralCoef < 1.5] <- "pink"
  color[df_spectralCoef >= 1.5] <- "brown"
  df_spectralCoef$color <- color
  
  group.colors <-
    c(
      blue = "#00008B",
      white = "white",
      pink = "#E30B5C",
      brown = "#A52A2A"
    )
  
  p_spectralCoef <-
    ggplot(df_spectralCoef, aes(x = beta, fill = color)) +
    geom_histogram(
      binwidth = binwidth,
      alpha = 0.5,
      position = "identity",
      color = 'black'
    ) +
    scale_fill_manual(values = group.colors) +
    
    geom_vline(
      aes(xintercept = 0),
      color = "black",
      linetype = "dashed",
      size = 1
    ) +
    geom_vline(
      aes(xintercept = 1),
      color = "black",
      linetype = "dashed",
      size = 1
    ) +
    geom_vline(
      aes(xintercept = 2),
      color = "black",
      linetype = "dashed",
      size = 1
    )
  p_spectralCoef  <-
    p_spectralCoef + theme(legend.position = "none") + theme_bw() + theme(
      plot.title = element_text(face = "bold", size = 15) ,
      axis.ticks = element_line(colour = "grey70", size = 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(size = 1, colour = 'grey70'),
      axis.line.y = element_line(size = 1, colour = 'grey70'),
      panel.border = element_blank(),
      legend.position = "bottom",
      text = element_text(size = 15)
    )
  
  return(p_spectralCoef)
  
}


# plots the changepoints in mean and variance
framework.plotCP <-
  function(dataAll,
           numPP,
           modelTitle,
           typePenalty,
           typeCP,
           freq = 1,
           xlab = "time",
           ylab,
           sampleList = c()) {
    # get the changePoints
    
    if (typePenalty == 'MBIC') {
      if (typeCP == 'individual') {
        cpSegMean <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_MBIC, "ind", "Mean"),
            changePoints_SegVals
          ))
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_MBIC, "ind", "Var"),
            changePoints_SegVals
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_MBIC,
              "ind",
              "VarMean"
            ),
            changePoints_SegVals
          ))
      } else{
        cpSegMean <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_MBIC, "All", "Mean"),
            changePoints_SegVals
          ))
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_MBIC, "All", "Var"),
            changePoints_SegVals
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_MBIC,
              "All",
              "VarMean"
            ),
            changePoints_SegVals
          ))
        
      }
      
    } else {
      if (typeCP == 'individual') {
        cpSegMean <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_AIC, "ind", "Mean"),
            changePoints_SegVals
          ))
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_AIC, "ind", "Var"),
            changePoints_SegVals
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_AIC,
              "ind",
              "VarMean"
            ),
            changePoints_SegVals
          ))
      } else{
        cpSegMean <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_AIC, "All", "Mean"),
            changePoints_SegVals
          ))
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_AIC, "All", "Var"),
            changePoints_SegVals
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_AIC,
              "All",
              "VarMean"
            ),
            changePoints_SegVals
          ))
      }
      
    }
    
    data <- dataAll %>%
      dplyr::select(data, ppID) %>%
      unnest(cols = c(data))
    data <- as.data.frame(data)
    data$cpSegMean <- cpSegMean
    data$cpSegVar <- cpSegVar
    data$cpSegVarMean <- cpSegVarMean
    
    if (length(sampleList) == 0) {
      results <-
        filter(data, ppID %in% sample(unique(data$ppID), numPP)) #select random participants for plotting
    } else{
      results <- filter(data, ppID %in% sampleList)
    }
    
    results$ppID <- as.factor(results$ppID)
    pModel <- ggplot(data = results)
    
    pModel <- pModel +
      geom_line(
        aes(x = time, y = y, group = ppID),
        color = "#616161",
        size = 0.4,
        alpha = 0.7
      )
    pModel <-
      pModel + labs(y = ylab, x = xlab) + scale_x_continuous(breaks =
                                                               seq(0, nObs, freq), labels = c(0:(nObs / freq)))
    # add segments for the mean
    pModel <-
      pModel + geom_line(aes(x = time, y = cpSegMean, group = ppID),
                         size = 0.8,
                         color = "firebrick")
    # add segments for the variance
    pModel <-
      pModel + geom_ribbon(
        aes(
          x = time,
          ymin = cpSegVarMean - cpSegVar * 2,
          ymax = cpSegVarMean + cpSegVar * 2,
          group = ppID
        ),
        fill = "grey70",
        alpha = 0.3
      )
    
    
    pModel <-
      pModel + labs(title = paste(modelTitle, "(" , typePenalty , ",", typeCP, ")")) #+ type + ")"
    pModel <-
      pModel +       theme(legend.position = "none") + theme_bw() + theme(
        plot.title = element_text(face = "bold", size = 15) ,
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey70', size = 1),
        legend.position = c(0.06, 0.75),
        text = element_text(size = 15)
      )
    
    
    return(pModel)
  }




# plots the changepoints in mean and variance
framework.plotPublication <-
  function(dataAll,
           numPP,
           modelTitle,
           typePenalty,
           typeCP,
           freq = 1,
           xlab = "time",
           ylab,
           sampleList = c(),regionFlag = FALSE) {
    # get the changePoints
    
    if (typePenalty == 'MBIC') {
      if (typeCP == 'individual') {
        cpSegMean <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_MBIC, "ind", "Mean"),
            changePoints_SegVals
          ))
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_MBIC, "ind", "Var"),
            changePoints_SegVals
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_MBIC,
              "ind",
              "VarMean"
            ),
            changePoints_SegVals
          ))
      } else{
        cpSegMean <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_MBIC, "All", "Mean"),
            changePoints_SegVals
          ))
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_MBIC, "All", "Var"),
            changePoints_SegVals
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_MBIC,
              "All",
              "VarMean"
            ),
            changePoints_SegVals
          ))
        
      }
      
    } else {
      if (typeCP == 'individual') {
        cpSegMean <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_AIC, "ind", "Mean"),
            changePoints_SegVals
          ))
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_AIC, "ind", "Var"),
            changePoints_SegVals
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_AIC,
              "ind",
              "VarMean"
            ),
            changePoints_SegVals
          ))
      } else{
        cpSegMean <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_AIC, "All", "Mean"),
            changePoints_SegVals
          ))
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_AIC, "All", "Var"),
            changePoints_SegVals
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_AIC,
              "All",
              "VarMean"
            ),
            changePoints_SegVals
          ))
      }
      
    }
    data <- dataAll %>%
      dplyr::select(data, ppID) %>%
      unnest(cols = c(data))
    data <- as.data.frame(data)
    data$cpSegMean <- cpSegMean
    data$cpSegVar <- cpSegVar
    data$cpSegVarMean <- cpSegVarMean
    
    otherStatsData <-
      cbind(
        dataAll$Statistics_Mean$ppID,
        dataAll$Statistics_Mean$m_sd,
        dataAll$Statistics_Mean$m_ac1,
        dataAll$Statistics_Stationary$stat_ac1,
        dataAll$Statistics_Mean$m_apprEntropy,
        dataAll$Statistics_Mean$m_spectralCoef
      )
    
    otherStatsData <- as.data.frame(otherStatsData)
    names(otherStatsData) <-
      c('ppID','sd', 'ac', 'stat_ac', 'entropy', 'color')
    
    
    
    if (length(sampleList) == 0) {
      results <-
        filter(data, ppID %in% sample(unique(data$ppID), numPP)) #select random participants for plotting
    } else{
      results <- filter(data, ppID %in% sampleList)
    }
    
    
    results$ppID <- as.factor(results$ppID)
    
    pModel <- ggplot(data = results)
    
    pModel <- pModel +
      geom_line(
        aes(x = time, y = y, group = ppID),
        color = "black", 
        size = 0.6,
        alpha = 0.8
      )
    breaks = seq(0, nObs, freq)
    labels = rep("", length(breaks))
    labels[1] <- 2006
    labels[nObs/freq +1] <- nObs/freq + 2006

    pModel <-
      pModel + labs(y = paste(ylab,"\n"), x = paste("\n",xlab)) + scale_x_continuous(breaks =
                                                              breaks, labels = labels)
    # add segments for the mean
    pModel <-
      pModel + geom_line(aes(x = time, y = cpSegMean, group = ppID),
                         size = 0.9,
                         color = "#1E8449")
 
    # add segments for the variance
    pModel <-
      pModel + geom_ribbon(
        aes(
          x = time,
          ymin = cpSegVarMean - cpSegVar * 2,
          ymax = cpSegVarMean + cpSegVar * 2,
          group = ppID
        ),
        fill = "grey70",
        alpha = 0.3
      )
    
    

    pModel <-
      pModel +       theme(legend.position = "none") + theme_bw() + theme(
        plot.title = element_text(face = "bold", size = 15) ,
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey70', size = 1),
        legend.position = c(0.06, 0.75),
        text = element_text(size = 15), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
      )
    
    
    if(regionFlag){
      region_names <- setNames(levels(data$region), levels(data$ppID))
      
      pNew_Sep <- pModel + facet_rep_wrap(~ ppID, ncol = 2,labeller=as_labeller(region_names)) + theme(
        strip.background = element_blank(),strip.text = element_textbox_highlight(size = 12, face = "bold"))
      
    }else{
      pNew_Sep <- pModel + facet_rep_wrap(~ ppID, ncol = 2) + theme(
        strip.background = element_blank(),strip.text = element_textbox_highlight(size = 12, face = "bold"))
      
    }
    

    # add annotations sd
    labelVec = c()
    for(idx in 1:length(sampleList)){
      labelVec[idx] = paste("ac:", round(otherStatsData[sampleList[idx], ]$ac, 3))
    }

    dat_text <- data.frame(
      label = labelVec,
      ppID   = sampleList,
      x     = rep(12, length(sampleList)),
      y     = rep(-2, length(sampleList))
    )

    pNew_Sep <- pNew_Sep + geom_text(data    = dat_text,
                                     mapping = aes(x = x, y = y, label = label),size =4,fontface = "bold.italic")





    # add annotations entropy

    labelVec = c()
    for(idx in 1:length(sampleList)){
      labelVec[idx] = paste("entropy:", round(otherStatsData[sampleList[idx], ]$entropy, 3))
    }

    dat_text <- data.frame(
      label = labelVec,
      ppID   = sampleList,
      x     = rep(65, length(sampleList)),
      y     = rep(-2, length(sampleList))
    )

    pNew_Sep <- pNew_Sep + geom_text(data    = dat_text,
                                     mapping = aes(x = x, y = y, label = label),size =4, fontface = "bold.italic")


    # add annotations color of noise
    labelVec = c()
    for(idx in 1:length(sampleList)){
      labelVec[idx] = paste("color of noise:", round(otherStatsData[sampleList[idx], ]$color, 3))
    }
    dat_text <- data.frame(
      label = labelVec,
      ppID   = sampleList,
      x     = rep(130, length(sampleList)),
      y     = rep(-2, length(sampleList))
    )

    pNew_Sep <- pNew_Sep + geom_text(data    = dat_text,
                                     mapping = aes(x = x, y = y, label = label), size =4, fontface = "bold.italic")


    return(pNew_Sep)
  }

# plots the changepoints in mean and variance
framework.plotPublication2 <-
  function(dataAll,
           numPP,
           modelTitle,
           typePenalty,
           typeCP,
           freq = 1,
           xlab = "time",
           ylab,
           sampleList = c(), nObs, type) {
    # get the changePoints
    
    if (typePenalty == 'MBIC') {
      if (typeCP == 'individual') {
        cpSegMean <-
          rep(unlist(map_dbl(dataAll$data, function(data) {
            mean(data$y)
          })), each = nObs)
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_MBIC, "ind", "Var"),
            changePoints_SegVals
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_MBIC,
              "ind",
              "VarMean"
            ),
            changePoints_SegVals
          ))
      } else{
        cpSegMean <-
          rep(unlist(map_dbl(dataAll$data, function(data) {
            mean(data$y)
          })), each = nObs)
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_MBIC, "All", "Var"),
            changePoints_SegVals
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_MBIC,
              "All",
              "VarMean"
            ),
            changePoints_SegVals
          ))
        
      }
      
    } else {
      if (typeCP == 'individual') {
        cpSegMean <-
          rep(unlist(map_dbl(dataAll$data, function(data) {
            mean(data$y)
          })), each = nObs)
        
        
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_AIC, "ind", "Var"),
            changePoints_SegVals
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_AIC,
              "ind",
              "VarMean"
            ),
            changePoints_SegVals
          ))
      } else{
        cpSegMean <-
          rep(unlist(map_dbl(dataAll$data, function(data) {
            mean(data$y)
          })), each = nObs)
        
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_AIC, "All", "Var"),
            changePoints_SegVals
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_AIC,
              "All",
              "VarMean"
            ),
            changePoints_SegVals
          ))
      }
      
    }
    data <- dataAll %>%
      dplyr::select(data, ppID) %>%
      unnest(cols = c(data))
    data <- as.data.frame(data)
    data$cpSegMean <- cpSegMean
    data$cpSegVar <- cpSegVar
    data$cpSegVarMean <- cpSegVarMean
    
    # get fitted model data 
    
    if (type == "mean") {
      # plot mean
      dataModel <- dataAll %>%
        dplyr::select(lm_valsMean, ppID) %>%
        unnest(cols = c(lm_valsMean))
        
      
      
    } else{
      # plot variance
      dataModel <- dataAll %>%
        dplyr::select(lm_valsVar, ppID) %>%
        unnest(cols = c(lm_valsVar)) 
        
    }
      
      dataModel <- as.data.frame(dataModel)
      
    
    
    # add trend to the data
    data$trend <- dataModel$.fitted
      
    if (length(sampleList) == 0) {
      results <-
        filter(data, ppID %in% sample(unique(data$ppID), numPP)) #select random participants for plotting
    } else{
      results <- filter(data, ppID %in% sampleList)
    }
    
    
    results$ppID <- as.factor(results$ppID)

    pModel <- ggplot(data = results)
    
    pModel <- pModel +
      geom_line(
        aes(x = time, y = y, group = ppID),
        color = "black",
        size = 0.7,
        alpha = 0.9
      )
    
    breaks = seq(0, nObs, freq)
    labels = rep("", length(breaks))
    labels[1] <- 2006
    labels[nObs/freq +1] <- nObs/freq + 2006
    
    pModel <-
      pModel + labs(y = paste(ylab,"\n"), x = paste("\n",xlab)) + scale_x_continuous(breaks =
                                                               breaks, labels = labels)
    # add segments for the mean
    pModel <-
      pModel + geom_line(aes(x = time, y = cpSegMean, group = ppID),
                         size = 0.8,
                         color = "#F1C40F", linetype = "dashed")
    # add segments for the variance
    pModel <-
      pModel + geom_ribbon(
        aes(
          x = time,
          ymin = cpSegVarMean - cpSegVar * 2,
          ymax = cpSegVarMean + cpSegVar * 2,
          group = ppID
        ),
        fill = "grey70",
        alpha = 0.3
      )
    
    # add trend line 
    pModel <- pModel + geom_line(aes(x = time, y = trend, group = ppID),
                                 size = 0.8,
                                 color = "#1E8449")
    
    pModel <-
      pModel +       theme(legend.position = "none") + theme_bw() + theme(
        plot.title = element_text(face = "bold", size = 15) ,
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey70', size = 1),
        legend.position = c(0.06, 0.75),
        text = element_text(size = 15), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
      )
    
    
    
    pNew_Sep <- pModel + facet_rep_wrap(~ ppID, ncol = 2) + theme(
      strip.background = element_blank(),strip.text = element_textbox_highlight(size = 12, face = "bold"))
    
    
    return(pNew_Sep)
  }


framework.plotPublication3 <-
  function(dataAll,
           numPP,
           modelTitle,
           typePenalty,
           typeCP,
           freq = 1,
           xlab = "time",
           ylab,
           sampleList = c(), nObs, type) {
    # get the changePoints
    
    if (typePenalty == 'MBIC') {
      if (typeCP == 'individual') {
        cpSegMean <-
          rep(unlist(map_dbl(dataAll$data, function(data) {
            mean(data$yVar)
          })), each = nObs)
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_MBIC, "ind", "Var"),
            changePoints_SegValsVar
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_MBIC,
              "ind",
              "VarMean"
            ),
            changePoints_SegValsVar
          ))
      } else{
        cpSegMean <-
          rep(unlist(map_dbl(dataAll$data, function(data) {
            mean(data$yVar)
          })), each = nObs)
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_MBIC, "All", "Var"),
            changePoints_SegValsVar
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_MBIC,
              "All",
              "VarMean"
            ),
            changePoints_SegValsVar
          ))
        
      }
      
    } else {
      if (typeCP == 'individual') {
        cpSegMean <-
          rep(unlist(map_dbl(dataAll$data, function(data) {
            mean(data$yVar)
          })), each = nObs)
        
        
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_AIC, "ind", "Var"),
            changePoints_SegValsVar
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_AIC,
              "ind",
              "VarMean"
            ),
            changePoints_SegValsVar
          ))
      } else{
        cpSegMean <-
          rep(unlist(map_dbl(dataAll$data, function(data) {
            mean(data$yVar)
          })), each = nObs)
        
        cpSegVar <-
          unlist(pmap(
            list(dataAll$data, dataAll$changePoints_AIC, "All", "Var"),
            changePoints_SegValsVar
          ))
        cpSegVarMean <-
          unlist(pmap(
            list(
              dataAll$data,
              dataAll$changePoints_AIC,
              "All",
              "VarMean"
            ),
            changePoints_SegValsVar
          ))
      }
      
    }
    data <- dataAll %>%
      dplyr::select(data, ppID) %>%
      unnest(cols = c(data))
    data <- as.data.frame(data)
    data$cpSegMean <- cpSegMean
    data$cpSegVar <- cpSegVar
    data$cpSegVarMean <- cpSegVarMean
    
    # get fitted model data 
    
    if (type == "mean") {
      # plot mean
      dataModel <- dataAll %>%
        dplyr::select(lm_valsMean, ppID) %>%
        unnest(cols = c(lm_valsMean))
      
      
      
    } else{
      # plot variance
      dataModel <- dataAll %>%
        dplyr::select(lm_valsVar, ppID) %>%
        unnest(cols = c(lm_valsVar)) 
      
    }
    
    dataModel <- as.data.frame(dataModel)
    
    
    
    # add trend to the data
    data$trend <- dataModel$.fitted
    
    if (length(sampleList) == 0) {
      results <-
        filter(data, ppID %in% sample(unique(data$ppID), numPP)) #select random participants for plotting
    } else{
      results <- filter(data, ppID %in% sampleList)
    }
    
    
    results$ppID <- as.factor(results$ppID)
    
    pModel <- ggplot(data = results)
    
    pModel <- pModel +
      geom_line(
        aes(x = time, y = yVar, group = ppID),
        color = "black",
        size = 0.7,
        alpha = 0.9
      )
    
    breaks = seq(0, nObs, freq)
    labels = rep("", length(breaks))
    labels[1] <- 2006
    labels[nObs/freq +1] <- nObs/freq + 2006
    
    pModel <-
      pModel + labs(y = paste(ylab, "\n"), x = paste("\n",xlab)) + scale_x_continuous(breaks =
                                                               breaks, labels = labels)
    # add segments for the mean
    pModel <-
      pModel + geom_line(aes(x = time, y = cpSegMean, group = ppID),
                         size = 0.8,
                         color = "#F1C40F", linetype = "dashed")
    # add segments for the variance
    pModel <-
      pModel + geom_ribbon(
        aes(
          x = time,
          ymin = cpSegVarMean - cpSegVar * 2,
          ymax = cpSegVarMean + cpSegVar * 2,
          group = ppID
        ),
        fill = "grey70",
        alpha = 0.3
      )
    
    # add trend line 
    pModel <- pModel + geom_line(aes(x = time, y = trend, group = ppID),
                                 size = 0.8,
                                 color = "#1E8449")
    
    pModel <-
      pModel +       theme(legend.position = "none") + theme_bw() + theme(
        plot.title = element_text(face = "bold", size = 15) ,
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey70', size = 1),
        legend.position = c(0.06, 0.75),
        text = element_text(size = 15), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
      )
    
    
    pNew_Sep <- pModel + facet_rep_wrap(~ ppID, ncol = 2) + theme(
      strip.background = element_blank(),strip.text = element_textbox_highlight(size = 12, face = "bold"))
    
    return(pNew_Sep)
  }


framework.plotMinMax <- function(time, y, data, title) {
  maxIDXes = which(y == max(y))
  minIDXes = which(y == min(y))
  
  if (length(maxIDXes) > 1) {
    maxIDX = sample(maxIDXes, 1)
  } else{
    maxIDX = maxIDXes
  }
  
  if (length(minIDXes) > 1) {
    minIDX = sample(minIDXes, 1)
  } else{
    minIDX = minIDXes
  }
  
  maxData <- data[[maxIDX]]$y
  minData <- data[[minIDX]]$y
  df <- as_tibble(data.frame(time, maxData, minData))
  dfL <- df %>%
    select(time, maxData, minData) %>%
    gather(type, value,-time)
  
  
  dfL$type <- recode(dfL$type, maxData = 'max', minData = 'min')
  
  p <-  dfL %>% ggplot(aes(x = time, y = value)) +
    geom_line() +
    facet_grid( ~ type) +
    theme(legend.position = "none") + theme_bw() + labs(title = title) +
    theme(
      plot.title = element_text(face = "bold", size = 15) ,
      axis.ticks = element_line(colour = "grey70", size = 0.2),
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = 'grey70', size = 1),
      legend.position = c(0.06, 0.75) ,
      text = element_text(size = 15)
    )
  
  return(list(
    'plot'  = p,
    'maxIDX' = maxIDX,
    'minIDX' = minIDX
  ))
}
