######################################################################################
## Read projections of climatic variables of interest for both scenarios
## RCP 4.5 and RCP 8.5, and every decade from 2010 to 2100.
## Modify the extent of the input rasters to match the default extent for Catalonia
######################################################################################

read.climatic.vars <- function(work.path){
  
  library(raster)
  library(tidyverse)
  
  ## Mask of the study area
  load("inputlyrs/rdata/mask.rdata")

  ## Default extent of raster maps of Catalonia  
  extCat <- extent(c(250000, 540000, 4480000, 4760000))
  
  for(clim.scn in c("rcp45", "rcp85")){
    for(decade in seq(10,90,10)){
      
      cat(paste("Building: scenario", clim.scn, "- decade", decade), "\n")
      
      ## Read annual minimum temp, annual precip and annual solar radiation
      ## Change extend to match the default
      TEMP <- raster(paste0(work.path, "/inputlyrs/asc/", clim.scn, "/", decade, "/mnan.asc"))
      TEMP <- extend(TEMP, extCat)
      PRECIP <- raster(paste0(work.path, "/inputlyrs/asc/", clim.scn, "/", decade, "/plan.asc"))
      PRECIP <- extend(PRECIP, extCat)
      RAD <- raster(paste0(work.path, "/inputlyrs/asc/", clim.scn, "/", decade, "/Radan.asc"))
      RAD <- extend(RAD, extCat)
      
      ## Build a data frame with MASK, TEMP, PRECIP and RAD
      ## And keep only cells from CAT
      clim <- data.frame(cell.id=1:ncell(MASK), mask=MASK[], temp=TEMP[], precip=PRECIP[], rad=RAD[])
      clim  <-  clim[!is.na(clim$mask),]
      clim <- select(clim, cell.id, temp, precip, rad)
      save(clim, file=paste0("inputlyrs/rdata/climate_", clim.scn, "_", decade, ".rdata"))
      
    } 
  } 
  
}