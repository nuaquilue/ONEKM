######################################################################################3
## Build land.rdata with the initialization of the state variables
######################################################################################3


#################TO DO#####################
#INCLUDE OPTION TO MASK VARIABLES WITH BARCEONA AREA
###########################################

read.state.vars <- function(work.path){
  
  library(raster)
  
  cat("Reading initial state variables", "\n")

  ## Read initial state vars 
  ## When MEDFIRE will be coupled to IPM, maybe MEDFIRE do not need to track 
  ## BIOMASS and AGE state variables.
  LCF <- raster(paste0(work.path, "/Medfire/inputlyrs/asc/LCFspp_1km_31N-ETRS89.asc")) 
  BIOMASS <- raster(paste0(work.path, "/Medfire/inputlyrs/asc/Biomass_1km_31N-ETRS89.asc")) ##change to updated ones
  AGE <- raster("Medfire/inputlyrs/asc/ForestAge_1km_31N-ETRS89.asc") ##change to updated ones
  TSDIST <- raster(paste0(work.path, "/Medfire/inputlyrs/asc/tsdist1km_Oriol.asc")) ##change to updated ones
  
  ## Build data frame with
  ## 1. cell.id, 2. spp, 3. biomass, 4. age, 5. tsdist, 6. distype, and 7. tburnt
  land <- data.frame(cell.id=1:ncell(LCF), spp=LCF[], biom=BIOMASS[], age=AGE[], tsdist=TSDIST[])
  land <- land[!is.na(land$spp),]
  land$distype <- NA; land$distype[land$spp<=17] <- 0
  land$tburnt <- land$distype
  
  ## Save it
  save(land, file="Medfire/inputlyrs/rdata/land.rdata")
   
}