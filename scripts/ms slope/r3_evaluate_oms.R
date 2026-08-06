##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Find the lowest-cAIC models for each species and 
##  calcualte model predictions across the bs_grid
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Import data  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Hardcoded path of the output folder
om_dir <- "G:/My Drive/Bering_redesign/"
species_list <- read.csv(paste0(om_dir, "data/species_list.csv"))
bs_grid <- readRDS(file = paste0(om_dir, "data/data_processed/grid_bs_year.RDS"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Calculate models with the lowest cAIC for each species
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model_caic <- data.frame()

for (irow in 1:nrow(x = species_list)) { ## loop over species -- start 
  
  # extract species and region. Each species has a specific survey footprint
  ispp = species_list$SCIENTIFIC_NAME[irow]
  iregion = species_list$FOOTPRINT[irow]
  
  # directory of the species' model output
  idir <- paste0(om_dir, "om/", iregion, "/", ispp, "/") 
  
  # extract the different model types for a given species from the patterns
  # in the filenames in idir
  imodels <- dir(path = idir) |> 
    grep(pattern = "cAIC", value = TRUE) |> 
    sub(pattern = "[^_]*_(.*)\\.RDS", replacement = "\\1")
  
  for (imodel in imodels) { # loop over imodels -- start
    
    ## append the hessian check and the caic 
    model_caic <- rbind(
      model_caic,
      data.frame(
        species = ispp, 
        region = iregion, 
        model = imodel,
        hessian = readRDS(paste0(idir, "sanity_", imodel, ".RDS"))$hessian_ok,
        caic = readRDS(file = paste0(idir, "cAIC_", imodel, ".RDS") )
      ) 
    )
    
  } # loop over imodels -- end
} ## loop over species -- end

## calculate lowest-caic model for each species and save output
best_model <- model_caic |> 
  subset(subset = hessian == TRUE) |>
  split(f = ~ species) |>
  lapply(FUN = function(x) x[which.min(x$caic), ]) |>
  do.call(what = "rbind")

saveRDS(object = best_model, file = paste0(om_dir, "om/best_model.RDS"))
write.csv(x = best_model, 
          file = paste0(om_dir, "om/best_model.csv"), 
          row.names = FALSE )

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Predict across the interpolation grid
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d_gct <- array(dim = c(max(bs_grid$cell),
                       nrow(x = best_model),
                       length(1982:2022)), 
               dimnames = list(1:max(bs_grid$cell),
                               best_model$species,
                               1982:2022))

for (irow in 1:nrow(x = best_model)) {
  
  iregion <- best_model$region[irow]
  ispp <- best_model$species[irow]
  survey_footprint <- list(bs_all = c("EBS", "NBS", "BSS"),
                           bs_shelf = c("EBS", "NBS"),
                           bs_slope = c("BSS"))[[iregion]]
  
  temp_fit <- readRDS(paste0(om_dir, "om/", iregion, "/", ispp, 
                             "/fit_", best_model$model[irow], ".RDS" ))
  
  ## Years to include: depends on region and species
  years_include <- sort(x = unique(x = temp_fit$data$year))
  
  temp_grid <- subset(x = bs_grid,
                      subset = region %in% survey_footprint 
                      & year %in% years_include)
  
  ## Cold pool as a spatially varying coefficient, append to temp_grid
  cp_df <- coldpool::cold_pool_index
  cp_index <- cbind(cp_df,
                    env = scale(coldpool::cold_pool_index$AREA_LTE2_KM2)) |>
    transform(year = as.integer(x = YEAR)) |>
    subset(select = c(year, env))
  temp_grid <- merge(x = temp_grid, y = cp_index, by = "year")
  
  ## Log depth and scale, append to temp_grid 
  mean_logdepth <- mean(x = log(x = temp_fit$data$depth_m))
  sd_logdepth <- sd(x = log(x = temp_fit$data$depth_m))
  temp_grid$scaled_depth <- 
    (log(x = temp_grid$depth_m) - mean_logdepth) / sd_logdepth
  
  ## Scale roms-sbt, append to temp_grid
  mean_roms_sbt <- mean(x = temp_fit$data$roms_sbt_c)
  sd_roms_sbt <- sd(x = temp_fit$data$roms_sbt_c)
  temp_grid$scaled_roms_sbt <- 
    (temp_grid$roms_sbt_c - mean_roms_sbt) / sd_roms_sbt
  
  # predict to the grid and insert into d_gct
  temp_preds <- predict(temp_fit, newdata = temp_grid, type = "response")
  d_gct[cbind(temp_preds$cell, ispp, temp_preds$year)] <- temp_preds$est
  
  cat("Finished with", ispp, "\n")
}

## Save output
saveRDS(object = d_gct, file = paste0(om_dir, "om/d_gct.RDS"))
