######################################################################################
##
######################################################################################

forest.mgmt <- function(land, coord, clim, orography, t){
  
  ## Tracking
  cat("Forest Management")
  
  ## Read management rules
  mgmt.rules <- read.table("inputfiles/MgmtRules.txt", header=T)
  
  ## Eq. of basal area to volume
  eq.ba.vol <- read.table("inputfiles/EqBasalAreaVol.txt", header=T)
  # eq.ba.volbark <- read.table("inputfiles/EqBasalAreaVolWithBark.txt", header=T)
  
  ## Sawlog and Wood demands
  dmnd <- read.table(paste0("inputfiles/", file.dmnd.harvest, ".txt"), header=T)
  dmnd.sawlog <- dmnd$Sawlogs[t]
  dmnd.wood <- dmnd$Primary[t]
  
  ## Num of neighbours in a circular neighbourhood according to radius (radius is in pixels)
  ## Assume that the neighbourhood is a star, with the maximum number of pixels in the
  ## east-west or north-south direction is 2*radius + 1 (1 is the center cell).
  ## The num of pixels is sequentially: 3+1*2, 5+3*2+1*2, 7+5*2+3*2+1*2, ...
      # nneigh <- seq(3,41,2) + cumsum(seq(1,40,2)*2)
  
  ## Find cells usable for management: forest with slope.pctg <=30% and dist.path <= 500m
  suit.mgmt <- left_join(select(land, -tsdist, -distype), 
                         select(orography, cell.id, slope.pctg, dist.path), by="cell.id") %>%
                  filter(spp <= 13 & slope.pctg <= 30 & dist.path <= 500) %>%
                  filter(spp != 9) # exclude quercus suber, not managed for sawlogs neither wood
  

  ## Idem for "final" harvesting
  if(dmnd.sawlog>0){
    suit.harvest <- filter(suit.mgmt, spp %in% c(1:7,12:13)) %>%
      left_join(select(clim, cell.id, sqi), by="cell.id") %>% left_join(mgmt.rules, by = c("spp", "sqi")) %>%
      filter(biom/10>=minab.fin) %>% mutate(ba.extract=pctgextract.fin*biom/1000) %>%
      left_join(eq.ba.vol, by = "spp") %>% mutate(vol.extract=cx*ba.extract+cx2*ba.extract*ba.extract) %>%
      mutate(priority=(slope.pctg+1)*(dist.path+1)*(1/vol.extract))
    suit.harvest$vol.sawlog <- suit.harvest$vol.extract*runif(nrow(suit.harvest), 0.90, 0.95)
    suit.harvest$vol.wood <- suit.harvest$vol.extract - suit.harvest$vol.sawlog
    suit.harvest <- suit.harvest[order(suit.harvest$priority, decreasing=F),]  
    cum.vol <- cumsum(suit.harvest$vol.sawlog) 
    if(max(cum.vol)>=dmnd.sawlog)
      cell.id.harvest <- suit.harvest$cell.id[1:which(cum.vol>dmnd.sawlog)[1]]
    else
      cell.id.harvest <- suit.harvest$cell.id
    dmnd.sawlog <- dmnd.sawlog - sum(suit.harvest$vol.sawlog[suit.harvest$cell.id %in% cell.id.harvest])
    dmnd.wood <- dmnd.wood - sum(suit.harvest$vol.wood[suit.harvest$cell.id %in% cell.id.harvest])
    harvesting <- filter(suit.harvest, cell.id %in% cell.id.harvest) %>%
                  select(cell.id, spp, vol.sawlog, vol.wood) %>% mutate(sylvi=3)
  }
  
  ## Idem for "dissmenatory" harvesting
  if(dmnd.sawlog>0){
    suit.harvest <- filter(suit.mgmt, spp %in% c(1:7,12:13)) %>%
      left_join(select(clim, cell.id, sqi), by="cell.id") %>% left_join(mgmt.rules, by = c("spp", "sqi")) %>%
      filter(biom/10>=minab.diss) %>%
      mutate(ba.extract=pmin(biom/10-thab.diss, pctgextract.diss*biom/1000)) %>%
      left_join(eq.ba.vol, by = "spp") %>% mutate(vol.extract=cx*ba.extract+cx2*ba.extract*ba.extract) %>%
      mutate(priority=(slope.pctg+1)*(dist.path+1)*(1/vol.extract))
    suit.harvest$vol.sawlog <- suit.harvest$vol.extract*runif(nrow(suit.harvest), 0.65, 0.85)
    suit.harvest$vol.wood <- suit.harvest$vol.extract - suit.harvest$vol.sawlog
    suit.harvest <- suit.harvest[order(suit.harvest$priority, decreasing=F),]  
    cum.vol <- cumsum(suit.harvest$vol.sawlog) 
    if(max(cum.vol)>=dmnd.sawlog)
      cell.id.harvest <- suit.harvest$cell.id[1:which(cum.vol>dmnd.sawlog)[1]]
    else
      cell.id.harvest <- suit.harvest$cell.id
    dmnd.sawlog <- dmnd.sawlog - sum(suit.harvest$vol.sawlog[suit.harvest$cell.id %in% cell.id.harvest])
    dmnd.wood <- dmnd.wood - sum(suit.harvest$vol.wood[suit.harvest$cell.id %in% cell.id.harvest])
    harvesting <- rbind(harvesting,
                        filter(suit.harvest, cell.id %in% cell.id.harvest) %>%
                          select(cell.id, spp, vol.sawlog, vol.wood) %>% mutate(sylvi=2))
  }
  
  
  ## Find locations suitable for "preparatory" harvesting according to current biomass
  ## Prioritize locations according to slope, dist.path and volume
  ## Harvest as much sawlogs as needed to meet the demand
  ## Between 65% - 85% of harvesting goes for sawlogs, the remainder goes for wood
  ## In the prioritization we should include some spatial (neighbourhood) criteria 
  ## and try to cluster the interventions
  if(dmnd.sawlog>0){
    suit.harvest <- filter(suit.mgmt, spp %in% c(1:7,12:13)) %>%
                    left_join(select(clim, cell.id, sqi), by="cell.id") %>% left_join(mgmt.rules, by = c("spp", "sqi")) %>%
                    filter(biom/10>=minab.prep) %>%
                    mutate(ba.extract=pmin(biom/10-thab.prep, pctgextract.prep*biom/1000)) %>%
                    left_join(eq.ba.vol, by = "spp") %>% mutate(vol.extract=cx*ba.extract+cx2*ba.extract*ba.extract) %>%
                    mutate(priority=(slope.pctg+1)*(dist.path+1)*(1/vol.extract))
    suit.harvest$vol.sawlog <- suit.harvest$vol.extract*runif(nrow(suit.harvest), 0.65, 0.85)
    suit.harvest$vol.wood <- suit.harvest$vol.extract - suit.harvest$vol.sawlog
    suit.harvest <- suit.harvest[order(suit.harvest$priority, decreasing=F),]  
    cum.vol <- cumsum(suit.harvest$vol.sawlog) 
    if(max(cum.vol)>=dmnd.sawlog)
      cell.id.harvest <- suit.harvest$cell.id[1:which(cum.vol>dmnd.sawlog)[1]]
    else
      cell.id.harvest <- suit.harvest$cell.id
    dmnd.sawlog <- dmnd.sawlog - sum(suit.harvest$vol.sawlog[suit.harvest$cell.id %in% cell.id.harvest])
    dmnd.wood <- dmnd.wood - sum(suit.harvest$vol.wood[suit.harvest$cell.id %in% cell.id.harvest])
    harvesting <- rbind(harvesting,
                        filter(suit.harvest, cell.id %in% cell.id.harvest) %>%
                        select(cell.id, spp, vol.sawlog, vol.wood) %>% mutate(sylvi=1) )
  }
  
  
  
  ## Now harvest for Wood
  if(dmnd.wood>0){
    suit.harvest <- filter(suit.mgmt, spp %in% 8:10) %>%
                    left_join(select(clim, cell.id, sqi), by="cell.id") %>% left_join(mgmt.rules, by = c("spp", "sqi")) %>%
                    filter(biom/10>=minab.fin) %>% mutate(ba.extract=pctgextract.fin*biom/1000) %>%
                    left_join(eq.ba.vol, by = "spp") %>% mutate(vol.extract=cx*ba.extract+cx2*ba.extract*ba.extract) %>%
                    mutate(priority=(slope.pctg+1)*(dist.path+1)*(1/vol.extract))
    suit.harvest$vol.sawlog <- 0
    suit.harvest$vol.wood <- suit.harvest$vol.extract - suit.harvest$vol.sawlog
    suit.harvest <- suit.harvest[order(suit.harvest$priority, decreasing=F),]  
    cum.vol <- cumsum(suit.harvest$vol.extract) 
    if(max(cum.vol)>=dmnd.wood)
      cell.id.harvest <- suit.harvest$cell.id[1:which(cum.vol>dmnd.wood)[1]]
    else
      cell.id.harvest <- suit.harvest$cell.id
    dmnd.wood <- dmnd.wood - sum(suit.harvest$vol.extract[suit.harvest$cell.id %in% cell.id.harvest])
    harvesting <- rbind(harvesting,
                        filter(suit.harvest, cell.id %in% cell.id.harvest) %>%
                        select(cell.id, spp, vol.sawlog, vol.wood) %>% mutate(sylvi=4))
  }
  
  return(harvesting)
}
