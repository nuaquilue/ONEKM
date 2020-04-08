######################################################################################
##
######################################################################################

update.clim <- function(MASK, land, orography, decade, clim.scn, clim.mdl){
  
  library(tidyverse)

  ## Tracking
  cat("Update climate", "\n")

  ## Read coefficients of site quality
  site.quality.spp <- read.table("inputfiles/SiteQualitySpp.txt", header=T)
  site.quality.index <- read.table("inputfiles/SiteQualityIndex.txt", header=T)
  site.quality.shrub <- read.table("inputfiles/SiteQualityShrub.txt", header=T)
  
  ## Update temp and precip
  load(paste0("inputlyrs/rdata/sdm_", clim.scn, "_", clim.mdl, "_", decade, ".rdata"))
  load(paste0("inputlyrs/rdata/climate_", clim.scn, "_", clim.mdl, "_", decade, ".rdata"))
  
  ## Join land.cover.spp
  clim$spp <- land$spp
  clim$aspect <- orography$aspect
  clim$slope <- orography$slope
  
  ## Assign SDM according to current spp distribution 
  clim$sdm <- NA
  clim$sdm[clim$spp==1] <- sdm$sdm.phalepensis[clim$spp==1]
  clim$sdm[clim$spp==2] <- sdm$sdm.pnigra[clim$spp==2]
  clim$sdm[clim$spp==3] <- sdm$sdm.ppinea[clim$spp==3]
  clim$sdm[clim$spp==4] <- sdm$sdm.psylvestris[clim$spp==4]
  clim$sdm[clim$spp==5] <- sdm$sdm.ppinaster[clim$spp==5]
  clim$sdm[clim$spp==6] <- sdm$sdm.puncinata[clim$spp==6]
  clim$sdm[clim$spp==7] <- sdm$sdm.aalba[clim$spp==7]
  clim$sdm[clim$spp==8] <- sdm$sdm.qilex[clim$spp==8]
  clim$sdm[clim$spp==9] <- sdm$sdm.qsuber[clim$spp==9]
  clim$sdm[clim$spp==10] <- sdm$sdm.qfaginea[clim$spp==10]
  clim$sdm[clim$spp==11] <- sdm$sdm.qhumilis[clim$spp==11]
  clim$sdm[clim$spp==12] <- sdm$sdm.fsylvatica[clim$spp==12]
  clim$sdm[clim$spp==13] <- sdm$sdm.other[clim$spp==13]
  clim$sdm[clim$spp==14] <- 1  ## SDM of shrub is always 1

  ## Compute SQ and SQI
  clim <- select(clim, cell.id, spp, temp, precip, sdm, aspect, slope) %>% 
          left_join(site.quality.spp, by="spp") %>% left_join(site.quality.index, by="spp") %>% 
          mutate(aux=c0+c_mnan*temp+c2_mnan*temp*temp+c_plan*precip+c2_plan*precip*precip+c_aspect*ifelse(aspect!=1,0,1)+c_slope*slope/10) %>%
          mutate(sq=1/(1+exp(-1*aux))) %>% mutate(sqi=ifelse(sq<=p50, 1, ifelse(sq<=p90, 2, 3))) %>%
          select(cell.id, spp, temp, precip, sdm, sqi)
  ## SQI for shrubs
  sqi.shrub <- filter(clim, spp==14) %>% select(spp, temp, precip) %>% left_join(site.quality.shrub, by="spp") %>%
               mutate(aux.brolla=c0_brolla+c_temp_brolla*temp+c_temp2_brolla*temp*temp+c_precip_brolla*precip+c_precip2_brolla*precip*precip,
                      aux.maquia=c0_maquia+c_temp_maquia*temp+c_temp2_maquia*temp*temp+c_precip_maquia*precip+c_precip2_maquia*precip*precip,
                      aux.boix=c0_boix+c_temp_boix*temp+c_temp2_boix*temp*temp+c_precip_boix*precip+c_precip2_boix*precip*precip,
                      sq.brolla=1/(1+exp(-1*aux.brolla)), sq.maquia=1/(1+exp(-1*aux.maquia)), sq.boix=1/(1+exp(-1*aux.boix)),
                      sqest.brolla=scale(sq.brolla), sqest.maquia=scale(sq.maquia), sqest.boix=scale(sq.boix),
                      sqi=ifelse(sqest.brolla>=sqest.maquia & sqest.brolla>=sqest.boix, 1,
                           ifelse(sqest.maquia>=sqest.brolla & sqest.maquia>=sqest.boix, 2,
                             ifelse(sqest.boix>=sqest.brolla & sqest.boix>=sqest.maquia, 3, 0))) )
  clim$sqi[clim$spp==14] <- sqi.shrub$sqi
  
  return(clim=clim)
}