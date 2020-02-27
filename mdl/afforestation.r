######################################################################################
##
######################################################################################

afforestation <- function(land, coord, orography, clim, sdm){
  
  ## Tracking
  cat("Afforestation", "\n")
  
  ## Read species reproductive age and afforestation model
        # age.spp <- read.table("inputfiles/SppAges.txt", header=T)  --> all spp have mature age 30
  afforest.mdl <- unlist(read.table("inputfiles/AfforestMdl.txt", header=T))
  seed.pressure <- unlist(read.table("inputfiles/SppSeedPressure.txt", header=T))
  
  ## Read coefficients of site quality models
  site.quality.spp <- read.table("inputfiles/SiteQualitySpp.txt", header=T)
  site.quality.index <- read.table("inputfiles/SiteQualityIndex.txt", header=T)
  
  ## Num of neighbours in a circular neighbourhood according to radius (radius is in pixels)
  ## Assume that the neighbourhood is a star, with the maximum number of pixels in the
  ## east-west or north-south direction is 2*radius + 1 (1 is the center cell).
  ## The num of pixels is sequentially: 3+1*2, 5+3*2+1*2, 7+5*2+3*2+1*2, ...
  nneigh <- seq(3,41,2) + cumsum(seq(1,40,2)*2)
  
  ## Coordinates of killed cells and their closest neighbours (do not count for the cell itself)
  shrub.coord <- filter(land, tsdist>=20, spp==14) %>% select(cell.id) 
  shrub.coord <- data.frame(cell.id=shrub.coord[sample(1:nrow(shrub.coord), 10000, replace=F),1])
  shrub.coord <- left_join(shrub.coord, coord, by = "cell.id")
  neigh.id <- nn2(coord[,-1], shrub.coord[,-1],  searchtype="priority", k=nneigh[shrub.colon.rad]) 
  neigh.id <- neigh.id$nn.idx  # dim 10.000 x 61
  
  ## Count number of forest mature neigbhours within their climatic range
  neigh.spp <- matrix(land$spp[neigh.id[,-1]], nrow=nrow(neigh.id), ncol=ncol(neigh.id)-1) 
  for(i in 1:13){
    neigh.sdm <- matrix(sdm[neigh.id[,-1],i+1], nrow=nrow(neigh.id), ncol=ncol(neigh.id)-1) 
    neigh.spp[neigh.spp==i] <- neigh.spp[neigh.spp==i] * neigh.sdm[neigh.spp==i]
  }
  neigh.spp <- neigh.spp * matrix(land$age[neigh.id[,-1]]>=30, nrow=nrow(neigh.id), ncol=ncol(neigh.id)-1) 
  rm(neigh.sdm); gc()
  
  ## Apply the afforestation model
  oldneigh <- apply(neigh.spp, 1, count.forest)
  z = afforest.mdl[1] + afforest.mdl[2]*clim$rad[neigh.id[,1]] + afforest.mdl[3]*(clim$rad[neigh.id[,1]])^2 +
      afforest.mdl[4]*20 + afforest.mdl[5]*(20^2) +  afforest.mdl[6]*(orography$slope[neigh.id[,1]]/10) + 
      afforest.mdl[7]*((orography$slope[neigh.id[,1]]/10)^2) + afforest.mdl[8]*(oldneigh/(ncol(neigh.id)-1)) +
      afforest.mdl[9]*((oldneigh/(ncol(neigh.id)-1))^2)
  p = 1/(1+exp(-1*z)); p = 1-(1-p)^(1/16); rm(z)
  
  ## 
  nneigh <- t(apply(neigh.spp, 1, count.spp.narm)) * 
            matrix(c(sdm[neigh.id[,1],1+1], sdm[neigh.id[,1],1+2], sdm[neigh.id[,1],1+3], sdm[neigh.id[,1],1+4],
                   sdm[neigh.id[,1],1+5], sdm[neigh.id[,1],1+6], sdm[neigh.id[,1],1+7], sdm[neigh.id[,1],1+8],
                   sdm[neigh.id[,1],1+9], sdm[neigh.id[,1],1+10], sdm[neigh.id[,1],1+11], sdm[neigh.id[,1],1+12],
                   sdm[neigh.id[,1],1+13]), nrow=nrow(neigh.id), ncol=13) *
            matrix(seed.pressure, nrow=nrow(neigh.id), ncol=13, byrow=T)
  available.neigh <- apply(nneigh, 1, exist.forest)
  
  ## Select species if available
  x <- runif(length(p),0,1)<=p & available.neigh 
  new.spp <- data.frame(cell.id=coord$cell.id[neigh.id[,1]*x],
                        spp=apply(nneigh[x,], 1, select.spp) )
  
  ## Join climatic and orographic variables to compute sq and then sqi
  new.spp <- left_join(new.spp, select(clim, cell.id, temp, precip), by = "cell.id") %>% 
             left_join(select(orography, cell.id, aspect, slope), by = "cell.id") %>%
             left_join(site.quality.spp, by = "spp") %>% left_join(site.quality.index, by = "spp") %>% 
             mutate(aux=c0+c_mnan*temp+c2_mnan*temp*temp+c_plan*precip+c2_plan*precip*precip+c_aspect*ifelse(aspect!=1,0,1)+c_slope*slope/10) %>%
             mutate(sq=1/(1+exp(-1*aux))) %>% mutate(sqi=ifelse(sq<=p50, 1, ifelse(sq<=p90, 2, 3))) %>%
             select(cell.id, spp, sqi) %>% mutate(biom=0, age=1, sdm=1)
  
  return(new.spp)
   
}

