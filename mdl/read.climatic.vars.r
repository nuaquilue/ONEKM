######################################################################################
## Read projections of climatic variables of interest for both scenarios
## RCP 4.5 and RCP 8.5, and every decade from 2010 to 2100.
## Modify the extent of the input rasters to match the default extent for Catalonia
######################################################################################

#################TO DO#####################
#CREATE MASK OF BARCELONA
#LOOk FOR DEFAULT EXTENT OF BARCELONA
###########################################

adapt.climatic.vars <- function(work.path){
  
  ## Mask of the new study area
  load("inputlyrs/rdata/mask.rdata")
  #MASK <- MASK_BCN
  new.mask <- data.frame(cell.id=1:ncell(MASK), x=MASK[])
  
  for(clim.scn in c("rcp45", "rcp85")){
    for(clim.mdl in c("KNMI-RACMO22E_ICHEC-EC-EARTH",
                      "KNMI-RACMO22E_MOHC-HadGEM2-ES",
                      "SMHI-RCA4_CNRM-CERFACS-CNRM-CM5",
                      "SMHI-RCA4_MPI-M-MPI-ESM-LR",
                      "SMHI-RCA4_MOHC-HadGEM2-ES")){
      for(decade in seq(10,90,10)){ 
        
        print(paste("Building: scenario", clim.scn, "model", clim.mdl, "- decade", decade))
        
        ## Load by default clima df and save them with other names
        load(paste0(work.path, "/inputlyrs/rdata/climate_", clim.scn, "_", clim.mdl, "_", decade, ".rdata"))
        save(clim, file=paste0(work.path, "/inputlyrs/rdata/climate_", clim.scn, "_", clim.mdl, "_", decade, "_CAT.rdata"))
        
        ## Filter data for the new study area
        clim <- left_join(clim, new.mask, by="cell.id") %>% filter(x==1) %>% select(-x)
        
        ## Save it
        save(clim, file=paste0(work.path, "/inputlyrs/rdata/climate_", clim.scn, "_", clim.mdl, "_", decade, ".rdata"))
        
      }
    }
  }
  
}

read.climatic.vars2 <- function(work.path, model){
  model <- model
  library(tidyverse)
  library(raster) 
  select <- dplyr::select
  ##change model
  #model <-  "SMHI-RCA4_MOHC-HadGEM2-ES"
  ##change study area
  CAT <- T
  BCN <- F  

  if (CAT){
  ## Mask of the study area
  load("inputlyrs/rdata/mask.rdata")
  ## Default extent of raster maps of Catalonia 
  extCat <- extent(c(250000, 540000, 4480000, 4760000))
  }  else if (BCN){
  ## Mask of the study area
  load("inputlyrs/rdata/mask_bcn.rdata")
  MASK = MASK_BCN
  ## Default extent of raster maps of Barcelona
  extCat <- extent(c(250000, 540000, 4480000, 4760000))} ##CHANGE TO DEFAULT EXTENT FOR BARCELONA
  
  for(clim.scn in c("rcp45", "rcp85")){
    for(decade in seq(00,90,10)){
    	cat(paste("Building: scenario", clim.scn, "- decade", decade), "\n")
    	clim <- data.frame(cell.id=1:ncell(MASK), mask=MASK[])
    	for(var in c("PRCPTOT", "TNMM")){        
	        ## Read annual minimum temp, annual precip and annual solar radiation
	        ## Change extend to match the default
	        if(var == "TNMM"){ #PRCPTOT_rcp45_KNMI-RACMO22E_ICHEC-EC-EARTH_proj00_1000m
	          if(decade == 0){
	            TEMP <- raster(paste0(work.path, "/inputlyrs/asc/clima_to_oriol/", var,"_",  clim.scn, "_", model,"_proj00_1000m.asc"))
	          }
	          else{
	            TEMP <- raster(paste0(work.path, "/inputlyrs/asc/clima_to_oriol/", var,"_",  clim.scn, "_", model,"_proj", decade, "_1000m.asc"))
	            }
	          TEMP <- extend(TEMP, extCat)
	          clim$temp = TEMP[]
	          }
	        else{
	          if(decade == 0){
	            PRECIP <- raster(paste0(work.path, "/inputlyrs/asc/clima_to_oriol/", var,"_",  clim.scn, "_", model,"_proj00_1000m.asc"))
	          }
	          else{
	            PRECIP <- raster(paste0(work.path, "/inputlyrs/asc/clima_to_oriol/", var,"_",  clim.scn, "_", model,"_proj", decade, "_1000m.asc"))
	            }
	          PRECIP <- extend(PRECIP, extCat)
	          clim$precip = PRECIP[]
	        }
      }
  	## Build a data frame with MASK, TEMP and PRECIP
    ## And keep only cells from CAT
    #clim <- data.frame(cell.id=1:ncell(MASK), mask=MASK[], temp=TEMP[], precip=PRECIP[], rad=RAD[])
    clim  <-  clim[!is.na(clim$mask),]
    clim <- select(clim, cell.id, temp, precip)
    save(clim, file=paste0("inputlyrs/rdata/climate_", clim.scn, "_", decade, ".rdata"))
    } 
  } 
  
}
read.climatic.vars <- function(work.path){
  
  library(tidyverse)
  library(raster) 
  select <- dplyr::select
  ##change study area
  CAT <- T
  BCN <- F  

  if (CAT){
	  ## Mask of the study area
	  load("inputlyrs/rdata/mask.rdata")
	  ## Default extent of raster maps of Catalonia 
	  extCat <- extent(c(250000, 540000, 4480000, 4760000))}
  else{ if (BCN){
	  ## Mask of the study area
	  load("inputlyrs/rdata/mask_bcn.rdata")
	  MASK = MASK_BCN
	  ## Default extent of raster maps of Barcelona
	  extCat <- extent(c(250000, 540000, 4480000, 4760000))}} ##CHANGE TO DEFAULT EXTENT FOR BARCELONA
  
  for(clim.scn in c("rcp45", "rcp85")){
    for(decade in seq(10,90,10)){
      
      cat(paste("Building: scenario", clim.scn, "- decade", decade), "\n")
      
      ## Read annual minimum temp, annual precip and annual solar radiation
      ## Change extend to match the default
      TEMP <- raster(paste0(work.path, "/inputlyrs/asc/", clim.scn, "/", decade, "/mnan.asc"))
      TEMP <- extend(TEMP, extCat)
      PRECIP <- raster(paste0(work.path, "/inputlyrs/asc/", clim.scn, "/", decade, "/plan.asc"))
      PRECIP <- extend(PRECIP, extCat)
      
      ## Build a data frame with MASK, TEMP, PRECIP and RAD
      ## And keep only cells from CAT
      clim <- data.frame(cell.id=1:ncell(MASK), mask=MASK[], temp=TEMP[], precip=PRECIP[], rad=RAD[])
      clim  <-  clim[!is.na(clim$mask),]
      clim <- select(clim, cell.id, temp, precip, rad)
      save(clim, file=paste0("inputlyrs/rdata/climate_", clim.scn, "_", decade, ".rdata"))
      
    } 
  } 
  
}