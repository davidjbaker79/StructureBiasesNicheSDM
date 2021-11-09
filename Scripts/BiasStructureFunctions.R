#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Functions for StructureBiasesNichesSDM analysis
#       By David Baker
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Spatial plots of variables and species occupancy patterns 
#'
#' This is just a simple function to make it easier to plot environment
#' variable, poc surfaces, and species occupancy patterns.
#'
#' @param dat Dataframe with X and Y coordinates and a variable to be plotted as
#' a fill variable
#' @param fill_var A variable to be mapped to \code{fill}.
#' @param out_path A path to save the figure.
#' @param points (logical) Can just plot points for occupancy patterns.
#' @param colPal A colour palette (viridis or scico).
#'
#' @import ggplot2
#' @import viridis
#' @import scico
#'
#' @return A spatial plot of simulated data.
#'
plot.variables <-
  function(dat,
           fill_var,
           out_path,
           points = FALSE,
           colPal = "viridis") {
    if (!points) {
      p1 <- ggplot() +
        geom_raster(aes_string(x = "X",
                               y = "Y",
                               fill = fill_var),
                    data = dat)
      if (colPal == "viridis")
        p1 <- p1 + scale_fill_viridis_c()
      if (colPal == "scico")
        p1 <- p1 + scale_fill_scico()
    } else {
      p1 <- ggplot() +
        geom_point(
          aes_string(x = "X",
                     y = "Y"),
          colour = "yellow",
          size = 0.05,
          data = dat[dat$pa == 1, ]
        )
    }
    
    p1 <- p1 + theme_nothing() %+replace%
      theme(panel.background = element_rect(fill = "black")) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    
    ggsave(
      filename = paste0(out_path, ".png"),
      plot = p1,
      width = 4,
      height = 4
    )
    p1
    
  }

#' Generate environmental variables with spatial autocorrelation
#'
#' @param n Gives the square dimensions of the model arena (n x n)
#' as the number of cells * 10
#' @param phi A vector defining the strength of autocorrelation. The lenght of
#' the vector determines the number of surfaces to return.
#'
#' @import mgcv
#' @import scales
#' @import Rfast
#'
#' @return Dataframe holding the generated environmental data.
#'
generate.enviro <- function(n = 10, phi = c(0.1, 0.1)) {
  # Set up a square lattice region
  simgrid <- expand.grid(1:n, 1:n)
  n <- nrow(simgrid)
  
  # Set up distance matrix
  distance <- as.matrix(Rfast::Dist(simgrid))
  
  # Generate random variable
  env_df <- simgrid[, 1:2] - 0.5
  names(env_df) <- c("X", "Y")
  for (i in 1:length(phi)) {
    env_df[[paste0("V", i)]] <-
      scales::rescale(mgcv::rmvn(1, rep(0, n), exp(-phi[i] * distance)))
  }
  
  # Visualize results
  gd <- rasterFromXYZ(env_df)
  gd <- disaggregate(gd, 5, method = "bilinear")
  plot(gd)
  
  # Data frame
  env_df <- data.frame(coordinates(gd), getValues(gd))
  names(env_df) <- c("X", "Y", names(gd))
  env_df$id <- 1:nrow(env_df)
  env_df <- env_df[, c("id", "X", "Y", paste0("V", 1:length(phi)))]
  
  # Return
  env_df
  
}

#' Virtual species generating function
#'
#' @param Var1 The fine-scale SAC surface (logistic response).
#' @param K The upper asymptote of curve for Var1.
#' @param B The inflection point of curve for Var1.
#' @param A The slope of the curve for Var1.
#' @param Var2 The coarse-scale SAC surface (Gaussian response).
#' @param k The upper limit of curve for Var2.
#' @param b The mean value of curve for Var2.
#' @param a The standard deviation of curve for Var2.
#' @param Var3 = intermediate-scale SAC surface (Gaussian response).
#'
#' @import scales
#'
#' @return Dataframe holding the species probability of occurrence data.
#'
species.occ <- function(Var1,
                        K = 0.8,
                        B = 0.6,
                        A = -0.05,
                        Var2,
                        k = 0.8,
                        b = 0.5,
                        a = 0.1,
                        Var3) {
  #--- Logistic response
  poc_1 <- K / (1 + exp((Var1 - B) / A))
  poc_g <-
    sapply(seq(0, 1, 0.01), function(x)
      K / (1 + exp((x - B) / A)))
  
  #--- Gaussian response
  poc_2 <- k * exp(-((Var2 - b) ^ 2 / (2 * a ^ 2)))
  poc_g <-
    sapply(seq(0, 1, 0.01), function(x)
      k * exp(-((x - b) ^ 2 / (2 * a ^ 2))))
  
  #--- Gaussian response
  poc_3 <- k * exp(-((Var3 - 0.5) ^ 2 / (2 * 0.25 ^ 2)))
  poc_g <-
    sapply(seq(0, 1, 0.01), function(x)
      k * exp(-((x - 0.5) ^ 2 / (2 * 0.25 ^ 2))))
  
  #--- Combine responses
  poc <- round(poc_1 * poc_2 * poc_3, 3)
  
  #--- Rescale (x - from[1])/diff(from) * diff(to) + to[1]
  poc <- scales::rescale(poc, to = c(0, 1))
  
  #---Convert to PA
  pa <- rbinom(n = length(poc),
               size = 1,
               prob = poc)
  
  #--- Join with XY coordinates
  pca_dat <- data.frame(poc = poc, pa = pa)
  
}

#' Sampling bias generating function
#'
#' This function defines a logistic curve (postive or )
#'
#' @param dat Data frame with X and Y coordinates and the third column variable
#' that SAC is generated in reference to.
#' @param K The upper asymptote of the SAC curve.
#' @param B The inflection point of the SAC curve.
#' @param A The slope of the SAC curve.
#'
#' @import scales
#'
#' @return Dataframe holding with added bias weights (wgt) column.
#'
samp.bias.Landcov <- function(dat,
                              K = 1,
                              B = 0.005,
                              A = 0.05) {
  dat$wgt <- K / (1 + exp((dat[, ncol(dat)] - B) / A))
  dat$wgt <- scales::rescale(dat$wgt)
  dat
  
}

#' Species virtual sampling
#'
#' @param sp_dat Data frame return by \code{species.occ}.
#' @param bias Vector of weights returned by \code{samp.bias.Landcov}.
#' @param nPres Numeric indicating the number of occurrence records to be sampled.
#' @param detProb Numeric indicating the probability of detecting a species on a
#' survey visit.
#' @param fpr Numeric indicating the rate at which species are miss identified.
#' @param listProb Numeric indicating the probability that a complete list is returned.
#'
#' @return Data frame with species det column added.
#'
species.sampling <- function(sp_dat,
                             bias = NULL,
                             nPres = 200,
                             detProb = 0.5,
                             fpr = 0,
                             listProb = 1) {
  # Join bias
  sp_dat <- cbind(sp_dat, bias)
  
  # Prevalence
  pr <- sp_dat[sp_dat$pa == 1, ]
  prev <- nrow(pr) / nrow(sp_dat)
  n_absen <- nPres / prev
  
  # Assign nPres detections based on biased surface
  pr <- pr[sample(1:nrow(pr), nPres, prob = pr$bias),]
  
  # Assign nPres detections based on biased surface
  ab <- sp_dat[sp_dat$pa == 0, ]
  ab <- ab[sample(1:nrow(ab), n_absen, prob = ab$bias),]
  
  # Combine data
  sp_dat_s <- rbind(pr, ab)
  
  # Assign detections based on detection probability
  sp_dat_s$det <- sp_dat_s$pa * rbinom(nrow(sp_dat_s), 1, detProb)
  
  # Adjust detections based on probability of  completing list
  sp_dat_s$det <- sp_dat_s$det * rbinom(nrow(sp_dat_s), 1, listProb)
  
  # Adjust detections based on probability of incorrect id
  det_id <- (1 - sp_dat_s$pa) * rbinom(nrow(sp_dat_s), 1, fpr)
  sp_dat_s$det <-  sp_dat_s$det + det_id
  
  # Return
  sp_dat_s
  
}

#' Run a complete experiment calling run_exp calling \code{run_exp_i}
#'
#' @param xp Dataframe containing settings for the experiment run.
#' @param modType String indicating type of model to be run.
#' @param modVar Vector giving names of model variables.
#' @param nReps Number of replications to run for each treatment.
#' @param biasMeth Type of bias correction method to apply.
#'
#' @return NULL
#'
run_exp <- function(xp,
                    modType = c("detNondet", "detOnly"),
                    biasCor = biasCor,
                    modVar = c("V1", "V2"),
                    nReps = 100,
                    biasMeth = c("None", "thin_maj", "thin_all", "cluster"))
{
  for (i in 1:nrow(xp)) {
    run_exp_i(
      spName = xp[i, "spName"],
      biasType = xp[i, "biasType"],
      detProb = xp[i, "detProb"],
      nPres = xp[i, "nPres"],
      falPosRat = xp[i, "falPosRat"],
      listInclusion = xp[i, "listInclusion"],
      effort = xp[i, "effort"],
      rep = nReps,
      modVar = modVar,
      modType = modType,
      biasCor = biasCor,
      biasMeth = biasMeth
    )
    
  }
  
  return(NULL)
  
}

#' Run each SDM experiment run within an experiment
#'
#' @param spName Name of species (one of "sp_com_gen" or "sp_com_spec").
#' @param biasType = Name of bias type.
#' @param detProb Numeric indicating the probability of detecting a species on a
#' survey visit.
#' @param nPres Numeric indicating the number of occurrence records to be sampled.
#' @param falPosRat Numeric indicating the rate at which species are miss
#' identified.
#' @param listInclusion Numeric indicating the probability that a complete list
#' is returned.
#' @param effort Numeric indicating the amount of visits to a site (time spent).
#' @param rep Number of replications to run for each treatment.
#' @param modVar Vector giving names of model variables.
#' @param biasCor Logical indicating whether bias correction should be applied.
#' @param modType String indicating type of model to be run.
#' @param biasMeth Type of bias correction method to apply.
#'
#' @return NULL
#'
run_exp_i <-
  function(spName,
           biasType,
           detProb,
           nPres,
           falPosRat,
           listInclusion,
           effort,
           rep = 100,
           modVar = c("V1", "V2"),
           biasCor,
           biasMeth,
           modType) {
    #- Load species data
    spPath <-
      paste0("Outputs/Species_data_files/", spName, ".Rdata")
    sp_in <- get(load(spPath))
    
    #- Load bias data
    if (biasType == "NULL") {
      landuse.bias <- rep(1, nrow(sp_in))
    } else {
      biasPath <- paste0("Outputs/Bias_data_files/", biasType, ".Rdata")
      landuse.bias <- get(load(biasPath))
      landuse.bias <- landuse.bias$wgt
    }
    
    #- Set up snowfall
    sfInit(parallel = TRUE, cpus = 8)
    sfExport(
      list = c(
        "sp_in",
        "landuse.bias",
        "biasType",
        "species.sampling",
        "thin_maj",
        "thin_all",
        "cluster_bkgrd",
        "detProb",
        "nPres",
        "falPosRat",
        "listInclusion",
        "effort",
        "modVar",
        "modType",
        "biasCor",
        "biasType",
        "spName",
        "biasMeth"
      )
    )
    
    sfLibrary(Metrics)
    sfLibrary(tidyverse)
    sfLibrary(gtools)
    sfLibrary(speedglm)
    sfLibrary(raster)
    sfLibrary(data.table)
    
    #- Run for n reps x = 1
    rep_all <- sfLapply(1:rep, function(x) {
      #- Sample distribution
      spSamp <- species.sampling(
        sp_dat = sp_in,
        bias = landuse.bias,
        nPres = nPres,
        detProb = detProb,
        fpr = falPosRat,
        listProb = listInclusion
      )
      
      #- Summary data
      n_sites <- length(unique(spSamp$id))
      n_visits <- nrow(spSamp)
      n_det <- sum(spSamp$det)
      n_occ <- sum(spSamp$pa)
      
      #- Model structure
      polVar <- sapply(modVar, function(i)
        paste0("poly(", i, ",2)"))
      modForm <-
        as.formula(paste0("det ~ ", paste(polVar, collapse = '+')))
      
      if (modType == "detNondet") {
        #- Select species
        groupVars <- lapply(c("id", "X", "Y", modVar), as.symbol)
        spDat <- spSamp %>%
          group_by(!!!groupVars) %>%
          summarise(det = max(det)) %>%
          as.data.frame()
        
        #- Implement spatial sub-sampling
        if (biasCor) {
          if (biasMeth == "thin_all")
            spDat <- thin_all(sp_in, spDat)
          if (biasMeth == "thin_maj")
            spDat <- thin_maj(sp_in, spDat)
        }
        
        #- Number of data points used to fit model
        n_mod_fit <- nrow(spDat)
        
        #- Fit model
        mod_1 <- speedglm(modForm,
                          data = spDat,
                          family = quasibinomial())
        
      }
      
      if (modType == "detOnly") {
        #print(x)
        #- Select species
        groupVars <- lapply(c('X', 'Y', modVar), as.symbol)
        spDat <- spSamp %>%
          group_by(!!!groupVars) %>%
          summarise(det = max(det)) %>%
          filter(det == 1) %>%
          as.data.frame()
        
        #- Assign weighting for background selection
        sp_in$wgt <- 1
        if (biasCor) {
          if (biasMeth == "cluster")
            sp_in <- cluster_bkgrd(sp_in, spDat)
          if (biasMeth == "real")
            sp_in$wgt <- landuse.bias
        }
        
        #- Select background sample
        bkgd <- sp_in[sample(
          1:nrow(sp_in),
          size = 10000,
          replace = TRUE,
          prob = sp_in$wgt
        ), ]
        bkgd <- dplyr::select(bkgd, !!!groupVars)
        bkgd$det <- 0
        
        #- Combine
        spBkDat <- rbind(spDat, bkgd)
        # ggplot(spBkDat) +
        # geom_point(aes(x = X, y = Y, colour = as.factor(det))) +
        #    theme_bw() + scale_colour_brewer(palette = "Set2")
        
        # Weight
        detN <- nrow(spBkDat[spBkDat$det == 1, ])
        bkgN <- nrow(spBkDat[spBkDat$det == 0, ])
        wgt <- c(rep(1, detN), rep(detN / bkgN, bkgN))
        
        # Number of data points used to fit model
        n_mod_fit <- nrow(spBkDat)
        
        # Weighted binomial glm
        mod_1 <- speedglm(modForm,
                          data = spBkDat,
                          family = quasibinomial(),
                          weights = wgt)
        
      }
      
      #-- Evaluate against 'real' distribution
      mod_pred <- cbind(sp_in,
                        fit = predict(mod_1,
                                      newdata = sp_in,
                                      type = "response"))
      #- AUC
      AUC <- Metrics::auc(mod_pred$pa, mod_pred$fit)
      #- Correlation coefficients
      corP <- cor(x = mod_pred$poc ,
                  y = mod_pred$fit,
                  method = "pearson")[[1]]
      corS <- cor(x = mod_pred$poc ,
                  y = mod_pred$fit,
                  method = "spearman")[[1]]
      #- Prediction errors
      RMSE <- rmse(mod_pred$fit, mod_pred$poc)
      MAE <- mae(mod_pred$fit, mod_pred$poc)
      #- WarI statistics
      sumPOC <- sum(mod_pred$poc)
      sumFit <- sum(mod_pred$fit)
      schrD <-
        1 - sum(abs((mod_pred$poc / sumPOC) - (mod_pred$fit / sumFit))) / 2
      warI <-
        1 - sum((sqrt(mod_pred$poc / sumPOC) - sqrt(mod_pred$fit / sumFit)) ^ 2) /
        2
      
      #-- Outputs x = 1
      sp_out <- data.frame(
        rep = x,
        modType = modType,
        biasCor = biasCor,
        biasMeth = biasMeth,
        spName = spName,
        detProb = detProb,
        nPres = nPres,
        fpr = falPosRat,
        listProb = listInclusion,
        rep_visit = effort,
        biasType = biasType,
        n_sites = n_sites,
        n_visits = n_visits,
        n_det = n_det,
        n_occ = n_occ,
        n_mod_fit = n_mod_fit,
        AUC = AUC,
        corP = corP,
        corS = corS,
        schrD = schrD,
        warI = warI,
        RMSE = RMSE,
        MAE = MAE
      )
      
    })
    dat_out <- do.call(rbind, rep_all)
    
    fnam <- paste(
      spName,
      biasType,
      modType,
      "BiasAdj",
      biasCor,
      biasMeth,
      nPres,
      sub("\\.", "", detProb),
      falPosRat,
      listInclusion,
      effort,
      sep = "_"
    )
    save(
      dat_out,
      file = paste0("Outputs/Model_results/", fnam, ".Rdata"),
      compress = "xz"
    )
    
    sfStop()
    
    return(NULL)
    
  }

#' Distance weighted background clustering
#'
#' Create raster then calculate distance to presences
#'
#' @param sp_in Dataframe produced by \code{species.sampling} from within
#' \code{run_exp_i}
#' @param spDat Data with species' detection.
#'
#' @import raster
#' @import dplyr
#'
#' @return Dataframe with spatial adjusted data.
#'
cluster_bkgrd <- function(sp_in, spDat) {
  spAll <- sp_in
  spAll <- left_join(spAll,
                     spDat[, c("X", "Y", "det")],
                     by = c("X", "Y"))
  r <- rasterFromXYZ(spAll[, c("X", "Y", "id", "det")])
  rbuf <- r[[2]]
  rbuf[rbuf == 0] <- NA
  rbuf <- distance(rbuf)
  rbuf <- data.frame(id = getValues(r[[1]]), dist = getValues(rbuf))
  sp_in <- left_join(sp_in, rbuf, by = c("id"))
  sp_in$wgt <-
    sapply(sp_in$dist, function(x)
      1 / (1 + exp((x - 4) / 1.5)))
  sp_in
  
}

#' Spatial thinning - thin all
#'
#' @param sp_in Dataframe produced by \code{species.sampling} from within
#' \code{run_exp_i}
#' @param spDat Data with species' detection.
#'
#' @import data.table
#' @import dplyr
#'
#' @return Dataframe with spatial adjusted data.
#'
thin_all <- function(sp_in, spDat) {
  spAll <- sp_in
  spAll <-
    left_join(spAll, spDat[, c("X", "Y", "det")], by = c("X", "Y"))
  
  #- Blocks
  blkSize <- 10
  nCells <- length(unique(spAll$id))
  blockLen <- sqrt(nCells) / blkSize
  blks <- matrix(NA, nrow =  sqrt(nCells), ncol =  sqrt(nCells))
  for (i in 1:blockLen) {
    st_i <- (blockLen * i) - (blockLen - 1)
    fin_i <- st_i + (blockLen - 1)
    blk <- rep(sort(rep(seq(st_i, fin_i, 1), blkSize)), blkSize)
    blks[, (i * blkSize - (blkSize - 1)):(i * blkSize)] <- blk
  }
  spAll$blk <- c(blks)
  
  #- Spatial weights
  wgt <- setDT(spAll)[, .(wgt = sum(!is.na(det))), blk]
  wdf <- merge(spAll, wgt, by = "blk")
  
  #- Get rid of zero because they have nothing to sample anyway
  wdf <- wdf[wdf$wgt > 0,]
  wdf$wgt <-  1 / wdf$wgt
  
  #- Thin data
  wdf <- wdf[wdf$id %in% unique(spDat$id),]
  wdf <- wdf[rbinom(1:nrow(wdf), size = 1, prob = wdf$wgt) == 1,]
  wdf <- wdf[, c("id", "X", "Y", "V1", "V2", "det")]
  
}

#' Spatial thinning - thin majority class
#'
#' @param sp_in Dataframe produced by \code{species.sampling} from within
#' \code{run_exp_i}
#' @param spDat Data with species' detection.
#'
#' @import data.table
#' @import dplyr
#'
#' @return Dataframe with spatial adjusted data.
#'
thin_maj <- function(sp_in, spDat) {
  #- Identify majority class
  ones <- length(spDat$det[spDat$det == 1])
  zeros <- length(spDat$det[spDat$det == 0])
  mcl <- 0
  if (ones >= zeros)
    mcl <- 1
  
  #- Create raster of sample data
  spAll <-  sp_in
  spAll <-
    merge(spAll, spDat[, c("X", "Y", "det")], by = c("X", "Y"), all.x = TRUE)
  
  #- Blocks
  blkSize <- 10
  nCells <- length(unique(spAll$id))
  blockLen <- sqrt(nCells) / blkSize
  blks <- matrix(NA, nrow =  sqrt(nCells), ncol =  sqrt(nCells))
  for (i in 1:blockLen) {
    st_i <- (blockLen * i) - (blockLen - 1)
    fin_i <- st_i + (blockLen - 1)
    blk <- rep(sort(rep(seq(st_i, fin_i, 1), blkSize)), blkSize)
    blks[, (i * blkSize - (blkSize - 1)):(i * blkSize)] <- blk
  }
  spAll$blk <- c(blks)
  
  #- Spatial weights
  wgt <- spAll
  wgt$det[wgt$det == 1 - mcl] <- NA
  wgt$det[wgt$det == 0 + mcl] <- 1
  wgt <- setDT(wgt)[, .(wgt = sum(!is.na(det))), blk]
  wdf <- merge(spAll, wgt, by = "blk")
  
  #- Get rid of zero because they have nothing to sample anyway
  wdf <- wdf[wdf$wgt > 0,]
  wdf$wgt <-  1 / wdf$wgt
  
  #- Thin data
  wdf <- wdf[wdf$id %in% unique(spDat$id),]
  wdf <- wdf[wdf$det == mcl, ]
  wdf <- wdf[rbinom(1:nrow(wdf), size = 1, prob = wdf$wgt) == 1,]
  wdf <- wdf[, c("id", "X", "Y", "V1", "V2", "det")]
  
  # Join thinned zeros with ones
  wdf <- rbind(spDat[spDat$det == 1 - mcl, ], wdf)
  
}
