######################################################################################
##
######################################################################################

post.fire <- function(land, coord, orography, clim, sdm){
  
  ## Tracking
  cat("Post-fire regeneration", "\n")
  
  ## Read matrix of secondary species according to species - sqi
  ## and fire response trait per species
  secondary.spp <- read.table("inputfiles/SecondarySpp.txt", header=T)
  response.trait <- read.table("inputfiles/FireResponseTrait.txt", header=T)
  spp.ages <- read.table("inputfiles/SppAges.txt", header=T)
  
  ## Read coefficients of site quality models
  site.quality.spp <- read.table("inputfiles/SiteQualitySpp.txt", header=T)
  site.quality.index <- read.table("inputfiles/SiteQualityIndex.txt", header=T)
  site.quality.shrub <- read.table("inputfiles/SiteQualityShrub.txt", header=T)
  
  ## Num of neighbours in a circular neighbourhood according to radius (radius is in pixels)
  ## Assume that the neighbourhood is a star, with the maximum number of pixels in the
  ## east-west or north-south direction is 2*radius + 1 (1 is the center cell).
  ## The num of pixels is sequentially: 3+1*2, 5+3*2+1*2, 7+5*2+3*2+1*2, ...
  nneigh <- seq(3,41,2) + cumsum(seq(1,40,2)*2)

  ## Coordinates of high-intensity burnt forest cells the current time step that
  ## - it doesn't regenerate per se (fire functional trait), or
  ## - it is out of its climatic range, or
  ## - it is younger than the regeneration age
  burnt.cells <- filter(land, spp<14 & tsdist==0, distype==hfire) %>% left_join(response.trait, by="spp") %>%
                 left_join(spp.ages, by="spp") %>% left_join(select(clim, cell.id, sdm, sqi), by="cell.id") %>%
                 filter(age<=regener | sdm==0 | trait==0) %>% left_join(coord, by = "cell.id") 
  
  ## Only continue if there's any cell with change of spp dominance
  if(nrow(burnt.cells)>0){
    ## Coordinates of their closest neighbours (do not count for the cell itself)
    neigh.id <- nn2(coord[,-1], select(burnt.cells,x,y),  searchtype="priority", k=nneigh[spp.distrib.rad])
    neigh.id <- neigh.id$nn.idx
    neigh.spp <- data.frame(cell.id=coord$cell.id[neigh.id[,1]],
                            matrix(land$spp[neigh.id[,-1]], nrow=nrow(neigh.id), ncol=ncol(neigh.id)-1) )
    
    ## Count number of neighbors per spp, assume that always there's a shrub cell in the neighbourhood
    neigh.spp <- data.frame(cell.id=coord$cell.id[neigh.id[,1]],
                            t(apply(neigh.spp[,-1], 1, count.spp))>=1 )
    neigh.spp$X14 <- T
    
    ## For those cells that a transition must be done:
    ## Look up sqi data and sencondary species  (according to dominant spp and sqi), 
    ## then add sdm of all tree species and finally
    ## add the number of forest spp in the neighbourhood
    burnt.cells <- left_join(burnt.cells, secondary.spp, by = c("spp", "sqi")) %>% left_join(sdm, by = "cell.id") %>% 
                   left_join(neigh.spp, by = "cell.id")
  
    ## Select spp among available
    new.cohort <- data.frame(cell.id=burnt.cells$cell.id,
                             spp=apply(select(burnt.cells, phalepensis:shrub) * 
                                         select(burnt.cells, sdm.phalepensis:sdm.shrub) * 
                                         select(burnt.cells, X1:X14), 1, select.cohort), 
                             biom=0, sdm=1, age=1)
    
    ## Join climatic and orographic variables to compute sq and then sqi
    new.cohort <- left_join(new.cohort, select(clim, cell.id, temp, precip), by = "cell.id") %>% 
      left_join(select(orography, cell.id, aspect, slope), by = "cell.id") %>%
      left_join(site.quality.spp, by = "spp") %>% left_join(site.quality.index, by = "spp") %>% 
      mutate(aux=c0+c_mnan*temp+c2_mnan*temp*temp+c_plan*precip+c2_plan*precip*precip+c_aspect*ifelse(aspect!=1,0,1)+c_slope*slope/10) %>%
      mutate(sq=1/(1+exp(-1*aux))) %>% mutate(sqi=ifelse(sq<=p50, 1, ifelse(sq<=p90, 2, 3))) %>%
      select(cell.id, spp, temp, precip, biom, age, sdm, sqi)
    sqi.shrub <- filter(new.cohort, spp==14) %>% select(spp, temp, precip) %>% left_join(site.quality.shrub, by = "spp") %>%
      mutate(aux.brolla=c0_brolla+c_temp_brolla*temp+c_temp2_brolla*temp*temp+c_precip_brolla*precip+c_precip2_brolla*precip*precip,
             aux.maquia=c0_maquia+c_temp_maquia*temp+c_temp2_maquia*temp*temp+c_precip_maquia*precip+c_precip2_maquia*precip*precip,
             aux.boix=c0_boix+c_temp_boix*temp+c_temp2_boix*temp*temp+c_precip_boix*precip+c_precip2_boix*precip*precip,
             sq.brolla=1/(1+exp(-1*aux.brolla)), sq.maquia=1/(1+exp(-1*aux.maquia)), sq.boix=1/(1+exp(-1*aux.boix)),
             sqi=ifelse(sq.brolla>=sq.maquia & sq.brolla>=sq.maquia, 1,
                        ifelse(sq.maquia>=sq.brolla & sq.maquia>=sq.boix, 2,
                               ifelse(sq.boix>=sq.brolla & sq.boix>=sq.maquia, 3, 0))) )
    new.cohort$sqi[new.cohort$spp==14] <- sqi.shrub$sqi
    
    return(select(new.cohort, -temp, -precip))  
  }
  
  else
    return(burnt.cells)
}

