######################################################################################
###  fire.regime()
###
######################################################################################

fire.regime <- function(land, coord, orography, pigni, swc, clim.sever, t, 
                        burnt.cells, fintensity, fire.ids, fire.id=0, annual.burnt=0){
                        
  options(warn=-1)
  cat(paste0("Fires in SWC: ", ifelse(swc==1, "Wind.", ifelse(swc==2, "Heat.", 
                               ifelse(swc==3, "Regular.", "Prescribed.")))))

  ## Function to select items not in a vector
  `%notin%` <- Negate(`%in%`)
  
  ## Read and load input data
  load("inputlyrs/rdata/pfst.pwind.rdata")
  clim.severity <- read.table(paste0("inputfiles/", file.clim.severity, ".txt"), header=T)
  pctg.hot.days <- read.table(paste0("inputfiles/", file.pctg.hot.days, ".txt"), header=T)
  prob.hot <- read.table("inputfiles/ProbHot.txt", header=T)
  prob.conv <- read.table("inputfiles/ProbConv.txt", header=T)
  aba.dist <- read.table("inputfiles/AnnualBurntAreaDist.txt", header=T)
  fs.dist <- read.table("inputfiles/FireSizeDist.txt", header=T)
  fire.supp <- read.table(paste0("inputfiles/", file.fire.suppression, ".txt"), header=T)
  spp.flammability <- read.table("inputfiles/SppSpreadRate.txt", header=T)
  fst.sprd.weight <- read.table(paste0("inputfiles/", file.sprd.weight, ".txt"), header=T)
  
  
  ## To be sure that non-burnable covers do not burn (water, rock, urban), nor agriculture land
  ## under prescribed burns
  if(swc<4)
    i <- land$spp<=17        # 2.938.560
  else
    i <- land$spp<=15        # 1.937.915
  subland <- land[i,]   
  suborography <- orography[i,]
  source.supp <- data.frame(cell.id=subland$cell.id, source.supp.sprd=F, source.supp.fuel=F)
  
  ## Reset TrackFires data frame each run and swc
  track.fire <- data.frame(year=NA, swc=NA, clim.sever=NA, fire.id=NA, fst=NA, 
                           wind=NA, atarget=NA, aburnt.highintens=NA, 
                           aburnt.lowintens=NA, asupp.sprd=NA, asupp.fuel=NA)
  
  
  ## Wind direction between 4 neigbours
  ## Wind direction is coded as 0-N, 90-E, 180-S, 270-W
  default.neigh <- data.frame(x=c(290L,-290L, -1L, 1L),
                                  windir=c(180,0,270, 90),
                                  dist=c(1000,1000,1000,1000))
  default.nneigh <- nrow(default.neigh)
  
  
  ## Find either fixed or stochastic annual target area for wildfires
  if(swc<4){ 
    ## Fixed
    if(sum(clim.severity[clim.severity$year==t,2:4])>0){  #
      is.aba.fix <- T
      area.target <- clim.severity[clim.severity$year==t, swc+1]
    }
    ## Find stochastic annual burnt area
    else{ 
      is.aba.fix <- F
      if(clim.sever==1 & swc<=2){  # decide if climatic severity is extrem for 'wind' or 'heat' swc
        pctg <- pctg.hot.days[pctg.hot.days$year==t, swc+1]
        prob.extrem <- 1/(1+exp(-(prob.hot$inter[swc] + prob.hot$slope[swc]*pctg)))
        if(runif(1,0,100) <= prob.extrem) # extreme
            clim.sever <- 2
      }
      area.target <- round(min(200000, max(10, 
                           rlnorm(1, aba.dist$meanlog[aba.dist$clim==clim.sever & aba.dist$swc==swc],
                                    aba.dist$sdlog[aba.dist$clim==clim.sever & aba.dist$swc==swc])))) 
    }  
  }
  ## Find annual target area for prescribed burns
  else{
    if(!is.na(pb.target.area))
      area.target <- pb.target.area
    else{
      accum.burnt.area[2:7] <- accum.burnt.area[1:6]
      accum.burnt.area[1] <- annual.burnt
      area.target <- pmax(0,pb.convenient.area*7-sum(accum.burnt.area))
    }
  }  
  cat(paste(" Annual target area:", area.target), "\n")
  
  
  ## Update prob.igni according to swc
  pigni <- mutate(pigni, psft=p*pfst.pwind[,ifelse(swc==1,1,2)+1]) %>%
           filter(cell.id %in% subland$cell.id)
  
  
  ## Pre-select the coordinates of old Mediterranean vegetation, i.e.
  ## Pinus halepensis, Pinus nigra, and Pinus pinea of age >=30 years.
  ## to compute probability of being a convective fire
  old.forest.coord <- filter(subland, spp<=3 & age>=30) %>% select(cell.id) %>% left_join(coord, by = "cell.id")

  
  ## Start burn until annual area target is not reached
  while(area.target>0){
    
    ## ID for each fire event
    fire.id <- fire.id+1
    
    ## Select an ignition point, to then decide the fire spread type, the fire suppression level,
    ## the wind direction and the target fire size according to clim and fire spread type
    ## What if selected igni has already been burnt?? How can I control it? pigni$psft==0 of burnt cells??
    igni.id <- sample(pigni$cell.id, 1, replace=F, pigni$psft)
    
    ## Assign the fire spread type 
    if(swc==1 | swc==3)
      fire.spread.type <- swc
    else if(swc==4)
      fire.spread.type <- 3
    else{
      neighs <- nn2(coord[,-1], filter(coord, cell.id==igni.id)[,-1], searchtype="standard", k=100)
      nneigh <- sum(neighs$nn.dists[,]<=500)  #sqrt(2*500^2)
      old.neighs <- nn2(old.forest.coord[,-1], filter(coord, cell.id==igni.id)[,-1], searchtype="standard", k=100)
      old.nneigh <- sum(old.neighs$nn.dists[,]<=500) #sqrt(2*500^2)
      z <- filter(prob.conv, clim==clim.sever)$inter + filter(prob.conv, clim==clim.sever)$slope*(old.nneigh/nneigh)*100
      1/(1+exp(-z))
      fire.spread.type <- ifelse(runif(1,0,1)<=1/(1+exp(-z)),2,3)
    }

    ## According to the fire spread type, look at the weights of each factor on spread rate
    rpb <- fst.sprd.weight[1,fire.spread.type+1]
    wwind <- fst.sprd.weight[2,fire.spread.type+1]
    wslope <- fst.sprd.weight[3,fire.spread.type+1]
    wflam <- fst.sprd.weight[4,fire.spread.type+1]
    waspc <- fst.sprd.weight[5,fire.spread.type+1]
    
    ## Assign the fire suppression levels
    sprd.th <- filter(fire.supp, clim==clim.sever, fst==fire.spread.type)$sprd.th
    fuel.th <- filter(fire.supp, clim==clim.sever, fst==fire.spread.type)$fuel.th
    
    ## Assign the main wind direction according to the fire spread type
    ## Wind directions: 0-N, 45-NE, 90-E, 135-SE, 180-S, 225-SW, 270-W, 315-NE
    if(fire.spread.type==1)  # N, NW or W according to map
      fire.wind <- sample(c(0,315,270), 1, replace=F, p=filter(pfst.pwind, cell.id==igni.id)[4:6])
    if(fire.spread.type==2)  # S 80%, SW 10%, SE 10%
      fire.wind <- sample(c(180,225,135), 1, replace=F, p=c(80,10,10))
    if(fire.spread.type==3)  # any at random
      fire.wind <- sample(seq(0,315,45), 1, replace=F)
    spp.flam <- filter(spp.flammability, fst==fire.spread.type) %>% select(-fst)
    
    ## Derive target fire size from a power-law according to clima and fire.spread.type 
    ## Or prescribed size from a log-normal
    if(swc<4){
      log.size <- seq(1.7, 5, 0.01)
      log.num <- filter(fs.dist, clim==clim.sever, fst==fire.spread.type)$intercept +
        filter(fs.dist, clim==clim.sever, fst==fire.spread.type)$slope * log.size
      fire.size.target <- sample(round(10^log.size), 1, replace=F, prob=10^log.num)
    }
    else
      fire.size.target <- max(1,min(round(rlnorm(1,pb.mean,pb.sd)),100))
    ## Bound fire.size.target to not exceed remaining area.target
    if(fire.size.target>area.target)
      fire.size.target <- area.target
    
    ## Initialize tracking variables
    ## Ignition always burnt, and it does in high intensity when no-PB
    fire.front <- igni.id
    aburnt.lowintens <- ifelse(swc==4, 100, 0)
    aburnt.highintens <- ifelse(swc==4, 0, 100)
    asupp.sprd <- asupp.fuel <- 0
    n.lowsprd <- n.lowfuel <- 0
    burnt.cells <- c(burnt.cells, igni.id)
    visit.cells <- c(burnt.cells, igni.id) # to account for burnt cells in previous SWC
    supp.sprd.cells <- supp.fuel.cells <- numeric()
    fintensity <- c(fintensity, 1)
    fire.ids <- c(fire.ids, fire.id)
      
    ## Start speading from active cells (i.e. the fire front)
    while((aburnt.lowintens+aburnt.highintens+asupp.sprd+asupp.fuel)<fire.size.target){
      
      ## Build a data frame with the theoretical 12 (=default.nneigh) neighbours of cells in fire.front, 
      ## and add the per definition wind direction and the distance.
      ## Filter cells thathave not been visited yet.
      neigh.id <- data.frame(cell.id=as.integer(rep(fire.front, each=default.nneigh)+
                                                rep(default.neigh$x, length(fire.front))),
                             source.id=rep(fire.front, each=default.nneigh),
                             dist=rep(default.neigh$dist,length(fire.front)),
                             windir=rep(default.neigh$windir,length(fire.front)) ) %>%
                  filter(cell.id %notin% visit.cells) %>% 
                  left_join(filter(source.supp, cell.id %in% fire.front), by=c("source.id" ="cell.id"))  # look if source cell has been suppressed
      
      
      ## Now find those neighbours that are currenty in Catalonia
      ## is_inCpp returns the position of neigh.id$cell.id in the 'land' data.frame (not the cell.id)!
      neigh.in.land <- is_inCpp(neigh.id$cell.id, subland$cell.id)
      i.land.in.neigh <- unique(neigh.in.land[which(neigh.in.land!=-1)])
      
      ## For all neighbours, compute fire intenstiy and flammability factors
      ## fire intenstiy and flam will be NA for non burnable covers
      neigh.land <- subland[i.land.in.neigh,] %>%
                    mutate(fuel=ifelse(spp %in% c(15,16,17), 0.5,
                                   ifelse(spp==14, 0.01638*biom,  # or 0.01638???
                                     ifelse(age<=7, 0.2,
                                       ifelse(biom<200, 0.4,
                                         ifelse(biom<480, 0.95, 0.6)))))) %>%
                    left_join(spp.flam, by="spp") %>% mutate(flam=wflam*flam)
      
      ## Now, add to i.land.in.neigh the indexes (positions) of fire.front cells (in case these are not already there)
      ## Further on, we'll need to know the elevation of the fire.front cells.
      i.land.in.neigh <- unique(c(i.land.in.neigh, is_inCpp(fire.front, subland$cell.id)) )
      
      ## Retrieve the orography variables for fire.front and neigbhour cells, 
      ## and already compute aspect factor
      neigh.orography <- suborography[i.land.in.neigh,] %>%
                         mutate(aspc=waspc*ifelse(aspect==1, 0.1, ifelse(aspect==3, 0.9, ifelse(aspect==4, 0.4, 0.3))))
      
      ## Get spread rate by:
      ## Joining to the neig.id data.frame the neigh.land and keep only burnable neighs 
      ## Joining to this df, the neigh.orography to get the elevation of the source cells
      ## Joining to this df, the neigh.orography to get the elevation of the neighbour cells
      ## Computing slope and wind factors
      sprd.rate <- left_join(neigh.land, neigh.id, by="cell.id") %>%
                   left_join(select(neigh.orography, cell.id, elev), by=c("source.id"="cell.id")) %>%
                   left_join(select(neigh.orography, cell.id, elev, aspc), by="cell.id")  %>% 
                   mutate(dif.elev = elev.y-elev.x, 
                          slope = wslope * pmax(pmin(dif.elev/dist,0.5),-0.5)+0.5, 
                          wind = wwind * (ifelse(abs(windir-fire.wind)>180, 
                                            360-abs(windir-fire.wind), abs(windir-fire.wind)))/180) %>% 
                   mutate(sr=slope+wind+flam+aspc, pb=1+rpb*log(sr*fuel)) %>%
                   group_by(cell.id) %>% 
                   summarize(fire.id=fire.id, spp=mean(spp), age=mean(age), fuel=max(fuel),
                             source.supp.sprd = any(source.supp.sprd[which(sr == max(sr))]),
                             source.supp.fuel = any(source.supp.fuel[which(sr == max(sr))]),
                             sr=max(sr), fintens=max(sr*fuel), pb=max(pb)) %>%
                   mutate(tosupp.sprd=(fintens<=sprd.th), tosupp.fuel=(spp<=14 & age<=fuel.th))
      
      ## Count how many cells could be suppressed
      n.lowsprd <- n.lowsprd + sum(sprd.rate$tosupp.sprd)
      n.lowfuel <- n.lowfuel + sum(sprd.rate$tosupp.fuel)
      
      ## Now compute actual burn state (T or F) according to pb and suppress:
      supp.sprd <- (sprd.rate$tosupp.sprd & n.lowsprd>=accum.supp) | sprd.rate$source.supp.sprd
      supp.fuel <- (sprd.rate$tosupp.fuel & n.lowfuel>=accum.supp) | sprd.rate$source.supp.fuel
      sprd.rate$burn <- runif(nrow(sprd.rate),0,pb.upper.th)<=sprd.rate$pb & sprd.rate$pb>pb.lower.th
      source.supp$source.supp.sprd[source.supp$cell.id %in% sprd.rate$cell.id[sprd.rate$tosupp.sprd]] <- T
      source.supp$source.supp.fuel[source.supp$cell.id %in% sprd.rate$cell.id[sprd.rate$tosupp.fuel]] <- T
      filter(source.supp, cell.id %in% sprd.rate$cell.id)
      # sprd.rate   
      
      ## Mark that all these neighs have been visited (before breaking in case no burn)
      visit.cells <- c(visit.cells, sprd.rate$cell.id)
      
      ## If at least there's a burn cell, continue, otherwise, stop
      if(!any(sprd.rate$burn))
        break
      
      ## Avoid fire overshooting at last iteration: Only burn cells with higher pb
      temp.burnt <- sprd.rate[sprd.rate$burn, c("cell.id", "pb")]
      if((aburnt.lowintens+aburnt.highintens+asupp.fuel+asupp.sprd+nrow(temp.burnt))>fire.size.target){
        max.burnt <- fire.size.target - (aburnt.lowintens+aburnt.highintens+asupp.fuel+asupp.sprd)
        temp.burnt <- temp.burnt[order(temp.burnt$pb, decreasing = TRUE),]
        def.burnt <- temp.burnt$cell.id[1:max.burnt]
        sprd.rate$burn <- (sprd.rate$cell.id %in% def.burnt)
      }
      
      ## Mark the burnt cells, the suppressed, and the fire intensity for burnt cells
      burnt.cells <- c(burnt.cells, sprd.rate$cell.id[sprd.rate$burn & !sprd.rate$tosupp.sprd & !sprd.rate$tosupp.fuel])
      supp.sprd.cells <- c(supp.sprd.cells, sprd.rate$cell.id[sprd.rate$tosupp.sprd & !sprd.rate$tosupp.fuel & !sprd.rate$burn])
      supp.fuel.cells <- c(supp.fuel.cells, sprd.rate$cell.id[sprd.rate$tosupp.fuel])
      fintensity <- c(fintensity, sprd.rate$fintens[sprd.rate$burn & !sprd.rate$tosupp.sprd & !sprd.rate$tosupp.fuel])
      fire.ids <- c(fire.ids, sprd.rate$fire.id[sprd.rate$burn & !sprd.rate$tosupp.sprd & !sprd.rate$tosupp.fuel])
      
      ## Select the new fire front
      exclude.th <- min(max(sprd.rate$sr)-0.005,   ## 'mad' -> median absolute deviation
                        rnorm(1,mean(sprd.rate$sr[sprd.rate$burn], na.rm=T)-mad(sprd.rate$sr[sprd.rate$burn], na.rm=T)/2,
                              mad(sprd.rate$sr[sprd.rate$burn], na.rm=T)))
      fire.front <- sprd.rate$cell.id[sprd.rate$burn & sprd.rate$sr>=exclude.th]
      
      ## Increase area burnt in either high or low intensity (Prescribed burns always burnt in low intensity)
      aburnt.lowintens <- aburnt.lowintens + sum(sprd.rate$burn & sprd.rate$fintens<=ifelse(swc<4,fire.intens.th,100))
      aburnt.highintens <- aburnt.highintens + sum(sprd.rate$burn & sprd.rate$fintens>ifelse(swc<4,fire.intens.th,100))
      asupp.sprd <- asupp.sprd + sum(sprd.rate$tosupp.sprd & !sprd.rate$tosupp.fuel & !sprd.rate$burn)
      asupp.fuel <- asupp.fuel + sum(sprd.rate$tosupp.fuel & !sprd.rate$burn)
      
      ## In the case, there are no cells in the fire front, stop trying to burn.
      ## This happens when no cells have burnt in the current spreading step
      if(length(fire.front)==0)
        break
      
    } # while 'fire'
    
    ## Write info about this fire
    track.fire <- rbind(track.fire, data.frame(year=t, swc, clim.sever, fire.id, fst=fire.spread.type, 
                                               wind=fire.wind, atarget=fire.size.target, aburnt.highintens, 
                                               aburnt.lowintens, asupp.sprd, asupp.fuel))
    # cat(paste("Fire:", fire.id, "- aTarget:", fire.size.target, "- aBurnt:", aburnt.lowintens+aburnt.highintens,
    #           "- aSupp:", asupp.sprd+asupp.fuel), "\n")
    
    ## Update annual burnt area
    area.target <- area.target - (aburnt.lowintens + aburnt.highintens + asupp.sprd + asupp.fuel)
    
  }  #while 'year'
  
  return(list(burnt.cells=burnt.cells, fintensity=fintensity,
              fire.ids=fire.ids, track.fire=track.fire[-1,]))
}