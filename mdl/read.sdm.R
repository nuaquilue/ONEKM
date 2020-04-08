######################################################################################
###  read.sdm()
###
### L'ordre d'espècies segueix sent:
### 1"Pinus sylvestris",2"Pinus uncinata",3"Pinus pinea",4"Pinus halepensis",5"Pinus nigra",
### 6"Pinus pinaster",7"Abies alba",8"Quercus faginea",9"Quercus ilex",10"Quercus suber",
### 11"Quercus humilis",12"Fagus sylvatica",13"OTrees."
######################################################################################

read.sdm <- function(work.path, set){
  
  library(sp)
  library(raster)
  library(tidyverse)
  
  ## Mask of the study area
  load(paste0(work.path, "/inputlyrs/rdata/mask.rdata"))
  
  ## Order species in sdm data.frames
  species <- c("phalepensis",	"pnigra",	"ppinea",	"psylvestris",	"ppinaster",	"puncinata",
               "aalba", "qilex", "qsuber",	"qfaginea",	"qhumilis",	"fsylvatica",	"other"  )
  order.spp.sdm <- c(4, 5, 3, 1, 6, 2, 7, 9, 10, 8, 11, 12, 13)
  
  ## 
  for(clim.scn in c("rcp45", "rcp85")){
    for(clim.mdl in c(#"KNMI-RACMO22E_ICHEC-EC-EARTH",
      #"KNMI-RACMO22E_MOHC-HadGEM2-ES",
      #                "SMHI-RCA4_CNRM-CERFACS-CNRM-CM5",
      #               "SMHI-RCA4_MPI-M-MPI-ESM-LR",
                      "SMHI-RCA4_MOHC-HadGEM2-ES")){
      for(decade in seq(10,90,10)){ 
        
        print(paste("Building: SDM", clim.scn, "model", clim.mdl, "- decade", decade))
        
        ## Load list of SDMs
        load(paste0(work.path, "/inputlyrs/rdata/SDM/SDM_", clim.scn,"_proj", decade, ".rdata"))
        ## Build a data frame with MASK and SDMs per spp
        sdm  <- data.frame(cell.id=1:ncell(MASK), mask=MASK[])        
        ## Change resolution from 1 km to 100m and merge data to the 'sdm' data frame
        for(i in order.spp.sdm){
          #sdm.proj[[i]] <- disaggregate(sdm.proj[[i]], fact=c(10,10))
          sdm <- cbind(sdm, data.frame(spp=sdm.proj[[i]][]))
        }
        names(sdm)[3:(length(order.spp.sdm)+2)] <- paste0("sdm.", species)
        ## Keep only cells in Catalonia, rename columns and add shrub sdm
        sdm <- filter(sdm, !is.na(mask)) %>% select(-mask) %>% mutate(sdm.shrub=1)
        ## Save it as a .rdata
        save(sdm, file=paste0(work.path, "/inputlyrs/rdata/sdm_", set, "_", clim.scn, "_", clim.mdl, "_", decade, ".rdata"))    
        
      }
    }
  }
        
}