######################################################################################
##  Remove any reference to SDMs
######################################################################################


update.clim <- function(MASK, land, orography, decade, clim.scn, psdm){
  
  library(tidyverse)

  ## Tracking
  cat("Updating climatic variables and site quality index", "\n")

  ## Read coefficients of site quality
  site.quality.spp <- read.table("inputfiles/SiteQualitySpp.txt", header=T)
  site.quality.index <- read.table("inputfiles/SiteQualityIndex.txt", header=T)
  site.quality.shrub <- read.table("inputfiles/SiteQualityShrub.txt", header=T)
  
  ## Update temp, precip and radiation
  load(paste0("inputlyrs/rdata/climate_", clim.scn, "_", decade, ".rdata"))
  
  ## Join land.cover.spp, aspect and slope data
  clim$spp <- land$spp
  clim$aspect <- orography$aspect
  clim$slope <- orography$slope

  ## Assign SDM according to current spp distribution 
      # load(paste0("inputlyrs/rdata/sdm_", psdm, "p_", clim.scn, "_", decade, ".rdata"))
      # clim$sdm <- NA
      # clim$sdm[clim$spp==1] <- sdm$sdm.phalepensis[clim$spp==1]
      # clim$sdm[clim$spp==2] <- sdm$sdm.pnigra[clim$spp==2]
      # clim$sdm[clim$spp==3] <- sdm$sdm.ppinea[clim$spp==3]
      # clim$sdm[clim$spp==4] <- sdm$sdm.psylvestris[clim$spp==4]
      # clim$sdm[clim$spp==5] <- sdm$sdm.ppinaster[clim$spp==5]
      # clim$sdm[clim$spp==6] <- sdm$sdm.puncinata[clim$spp==6]
      # clim$sdm[clim$spp==7] <- sdm$sdm.aalba[clim$spp==7]
      # clim$sdm[clim$spp==8] <- sdm$sdm.qilex[clim$spp==8]
      # clim$sdm[clim$spp==9] <- sdm$sdm.qsuber[clim$spp==9]
      # clim$sdm[clim$spp==10] <- sdm$sdm.qfaginea[clim$spp==10]
      # clim$sdm[clim$spp==11] <- sdm$sdm.qhumilis[clim$spp==11]
      # clim$sdm[clim$spp==12] <- sdm$sdm.fsylvatica[clim$spp==12]
      # clim$sdm[clim$spp==13] <- sdm$sdm.other[clim$spp==13]
      # clim$sdm[clim$spp==14] <- 1  ## SDM of shrub is always 1
  clim$sdm <- 1

  ## Compute SQ and SQI for forest species
  clim <- select(clim, cell.id, spp, temp, precip, rad, sdm, aspect, slope) %>% 
          left_join(site.quality.spp, by="spp") %>% left_join(site.quality.index, by="spp") %>% 
          mutate(aux=c0+c_mnan*temp+c2_mnan*temp*temp+c_plan*precip+c2_plan*precip*precip+c_aspect*ifelse(aspect!=1,0,1)+c_slope*slope/10) %>%
          mutate(sq=1/(1+exp(-1*aux))) %>% mutate(sqi=ifelse(sq<=p50, 1, ifelse(sq<=p90, 2, 3))) %>%
          select(cell.id, spp, temp, precip, rad, sdm, sqi)
  
  ## SQI for shrubs
  sqi.shrub <- filter(clim, spp==14) %>% select(spp, temp, precip, rad) %>% left_join(site.quality.shrub, by="spp") %>%
               mutate(aux.brolla=c0_brolla+c_temp_brolla*temp+c_temp2_brolla*temp*temp+c_precip_brolla*precip+c_precip2_brolla*precip*precip,
                      aux.maquia=c0_maquia+c_temp_maquia*temp+c_temp2_maquia*temp*temp+c_precip_maquia*precip+c_precip2_maquia*precip*precip,
                      aux.boix=c0_boix+c_temp_boix*temp+c_temp2_boix*temp*temp+c_precip_boix*precip+c_precip2_boix*precip*precip,
                      sq.brolla=1/(1+exp(-1*aux.brolla)), sq.maquia=1/(1+exp(-1*aux.maquia)), sq.boix=1/(1+exp(-1*aux.boix)),
                      sqi=ifelse(sq.brolla>=sq.maquia & sq.brolla>=sq.maquia, 1,
                           ifelse(sq.maquia>=sq.brolla & sq.maquia>=sq.boix, 2,
                             ifelse(sq.boix>=sq.brolla & sq.boix>=sq.maquia, 3, 0))) )
  clim$sqi[clim$spp==14] <- sqi.shrub$sqi
  
  
  return(clim=clim)
}