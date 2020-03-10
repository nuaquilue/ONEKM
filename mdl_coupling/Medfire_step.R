Medfire.step <- function(t){

  load("./mdl_coupling/output/land.rdata")
    
  ## Track scenario, replicate and time step
  cat(paste0("scn: ", scn.name," - run: ", irun, "/", nrun, " - time: ", t, "/", time.horizon), "\n")
  
  if(write.sp.outputs)
    BURNT <- MASK
  
  ## 1. CLIMATE CHANGE  
  if(processes[clim.id] & t %in% temp.clim.schedule){
    ## Update temp, precip, solar radiation and site quality index
    clim <- update.clim(MASK, land, orography, decade=(1+floor(t/10))*10, clim.scn)
    temp.clim.schedule <- temp.clim.schedule[-1] 
  }
  
  
  ## 2. LAND-COVER CHANGE
  chg.cells <- integer()
  if(processes[lchg.id] & t %in% temp.lchg.schedule){
    lchg.out <- land.cover.change(land, coord, orography, lc.trans=1, chg.cells)
    # update interface values
    interface <- update.interface(land)
    temp.lchg.schedule <- temp.lchg.schedule[-1] 
  }
  
  
  ## 3. FOREST MANAGEMENT
  if(processes[fmgmt.id] & t %in% temp.fmgmt.schedule){
    aux <- forest.mgmt(land, coord, clim, harvest, t)
    land$tsdist[land$cell.id %in% aux$cell.id] <- 0
    land$distype[land$cell.id %in% aux$cell.id] <- fmgmt.id*10+aux$sylvi
    track.fmgmt <- rbind(track.fmgmt, 
                         data.frame(run=irun, year=t, 
                                    group_by(aux, spp, sylvi) %>% summarize(sawlog=sum(vol.sawlog), wood=sum(vol.wood))))
    temp.fmgmt.schedule <- temp.fmgmt.schedule[-1] 
  }
  
  
  ## 4. FIRE
  ## Tracking variables to be re-initialized each time step
  ## Out of the "if(fires)" in case only prescribed burns are applied
  burnt.cells <- integer()
  burnt.intens <- integer()
  annual.burnt <- 0
  if(processes[fire.id] & t %in% temp.fire.schedule){
    pigni <- prob.igni(land, orography, clim, interface)
    # Decide climatic severity of the year (default is mild)
    clim.sever <- 0
    if(runif(1,0,100) < clim.severity[clim.severity$year==t, ncol(clim.severity)]) # not-mild
      clim.sever <- 1
    # swc = wind
    fire.out <- fire.regime(land, coord, orography, pigni, swc=1, clim.sever, t, 
                            burnt.cells, burnt.intens, annual.burnt)
    burnt.cells <- fire.out[[1]]; burnt.intens <- fire.out[[2]]
    if(nrow(fire.out[[3]])>0)
      track.fire <- rbind(track.fire, data.frame(run=irun, fire.out[[3]]))
    annual.burnt <- annual.burnt+sum(fire.out[[3]]$aburnt.highintens + fire.out[[3]]$aburnt.lowintens)
    # swc = heat
    fire.out <- fire.regime(land, coord, orography, pigni, swc=2, clim.sever, t, 
                            burnt.cells, burnt.intens, annual.burnt)
    burnt.cells <- fire.out[[1]]; burnt.intens <- fire.out[[2]]
    if(nrow(fire.out[[3]])>0)
      track.fire <- rbind(track.fire, data.frame(run=irun, fire.out[[3]]))
    annual.burnt <- annual.burnt+sum(fire.out[[3]]$aburnt.highintens + fire.out[[3]]$aburnt.lowintens)
    # swc = regular
    fire.out <- fire.regime(land, coord, orography, pigni, swc=3, clim.sever, t, 
                            burnt.cells, burnt.intens, annual.burnt)
    burnt.cells <- fire.out[[1]]; burnt.intens <- fire.out[[2]]
    if(nrow(fire.out[[3]])>0)
      track.fire <- rbind(track.fire, data.frame(run=irun, fire.out[[3]]))
    annual.burnt <- annual.burnt+sum(fire.out[[3]]$aburnt.highintens + fire.out[[3]]$aburnt.lowintens)
    # done with fires!
    land$tsdist[land$cell.id %in% burnt.cells] <- 0
    land$tburnt[land$cell.id %in% burnt.cells] <- land$tburnt[land$cell.id %in% burnt.cells] + 1
    land$distype[land$cell.id %in% burnt.cells[as.logical(burnt.intens)]] <- hfire
    land$distype[land$cell.id %in% burnt.cells[!as.logical(burnt.intens)]] <- lfire
    temp.fire.schedule <- temp.fire.schedule[-1] 
    rm(fire.out)
    save(burnt.cells, file=paste0(out.path, "/burnt_cells.rdata"))
  }
  
  
  ## 5. PRESCRIBED BURNS
  if(processes[pb.id] & t %in% temp.pb.schedule){
    fire.out <- fire.regime(land, coord, orography, pigni, swc=4, clim.sever, t, 
                            burnt.cells, burnt.intens, annual.burnt)
    pb.cells <- fire.out[[1]]
    if(nrow(fire.out[[3]])>0)
      track.pb <- rbind(track.pb, data.frame(run=irun, fire.out[[3]][,c(1,3,4,6,7,9)]))
    # done with prescribed burns!
    pb.cells <- pb.cells[!(pb.cells %in% burnt.cells)]
    land$tsdist[land$cell.id %in% pb.cells] <- 0
    land$tburnt[land$cell.id %in% pb.cells] <- land$tburnt[land$cell.id %in% pb.cells] + 1
    land$distype[land$cell.id %in% pb.cells] <- pb
    temp.pb.schedule <- temp.pb.schedule[-1] 
    rm(fire.out)
  }
  
  
  ## 6. DROUGHT
  killed.cells <- integer()
  if(processes[drought.id] & t %in% temp.drought.schedule){
    killed.cells <- drought(land, clim, t)
    land$tsdist[land$cell.id %in% killed.cells] <- 0
    land$distype[land$cell.id %in% killed.cells] <- drght
    track.drougth <- rbind(track.drougth,
                          data.frame(run=irun, year=t, 
                                     filter(land, cell.id %in% killed.cells) %>% group_by(spp) %>% summarize(ha=length(spp))) )
    temp.drought.schedule <- temp.drought.schedule[-1] 
  }
  
  
  ## 7. POST-FIRE REGENERATION
  if(processes[post.fire.id] & t %in% temp.post.fire.schedule){
    aux  <- post.fire(land, coord, orography, clim, sdm)
    if(nrow(aux)>0){
      spp.out <- land$spp[land$cell.id %in% aux$cell.id]
      land$spp[land$cell.id %in% aux$cell.id] <- aux$spp
      land$biom[land$cell.id %in% aux$cell.id] <- growth.10y(aux, aux)
      clim$spp[clim$cell.id %in% aux$cell.id] <- aux$spp
      clim$sdm[clim$cell.id %in% aux$cell.id] <- 1
      clim$sqi[clim$cell.id %in% aux$cell.id] <- aux$sqi
      track.post.fire <- rbind(track.post.fire, data.frame(run=irun, year=t, table(spp.out, aux$spp)))  
    }
    rm(aux) #; rm(burnt.cells); rm(pb.cells)
    temp.post.fire.schedule <- temp.post.fire.schedule[-1] 
  }
  
  
  ## 8. COHORT ESTABLISHMENT
  if(processes[cohort.id] & t %in% temp.cohort.schedule & length(killed.cells)>0){
    aux  <- cohort.establish(land, coord, orography, clim, sdm)
    spp.out <- land$spp[land$cell.id %in% killed.cells]
    land$spp[land$cell.id %in% killed.cells] <- aux$spp
    land$biom[land$cell.id %in% killed.cells] <- growth.10y(aux, aux)
    clim$spp[clim$cell.id %in% killed.cells] <- aux$spp
    clim$sdm[clim$cell.id %in% killed.cells] <- 1
    clim$sqi[clim$cell.id %in% killed.cells] <- aux$sqi
    track.cohort <- rbind(track.cohort, data.frame(run=irun, year=t, table(spp.out, aux$spp)))
    rm(aux); rm(killed.cells); gc()
    temp.cohort.schedule <- temp.cohort.schedule[-1] 
  }
  
  
  ## 9. AFFORESTATION
  if(processes[afforest.id] & t %in% temp.afforest.schedule){
    aux  <- afforestation(land, coord, orography, clim, sdm)
    land$spp[land$cell.id %in% aux$cell.id] <- aux$spp
    land$biom[land$cell.id %in% aux$cell.id] <- growth.10y(aux, aux)
    clim$spp[clim$cell.id %in% aux$cell.id] <- aux$spp
    clim$sdm[clim$cell.id %in% aux$cell.id] <- 1
    clim$sqi[clim$cell.id %in% aux$cell.id] <- aux$sqi
    track.afforest <- rbind(track.afforest, data.frame(run=irun, year=t, table(aux$spp)))
    temp.afforest.schedule <- temp.afforest.schedule[-1] 
  }
  
  
  ## 10. GROWTH
  if(processes[growth.id] & t %in% temp.growth.schedule){
    land$biom <- growth.10y(land, clim)
    land$age <- pmin(land$age+1,600)
    land$tsdist <- pmin(land$tsdist+1,600)
    aux <- filter(land, spp<=13) %>% select(spp, biom) %>% left_join(eq.ba.vol, by="spp") %>% 
           mutate(vol=cx*biom/10+cx2*biom*biom/100) %>% select(-cx, -cx2) %>%
           left_join(eq.ba.volbark, by="spp") %>% 
           mutate(volbark=cx*biom/10+cx2*biom*biom/100) %>% select(-cx, -cx2) %>% 
           left_join(eq.ba.carbon, by="spp") %>% 
           mutate(carbon=c*biom/10) %>% group_by(spp) %>% select(-c) %>%
           summarise(area=length(vol), vol=sum(vol), volbark=sum(volbark), carbon=sum(carbon))  
    aux.shrub <- filter(land, spp==14) %>% select(spp, biom) %>% group_by(spp) %>%
                 summarise(area=length(biom), vol=sum(biom), volbark=0, carbon=0)  
    track.land <- rbind(track.land, data.frame(run=irun, year=t, aux), data.frame(run=irun, year=t, aux.shrub))
    temp.growth.schedule <- temp.growth.schedule[-1] 
  }
  
  ## Print maps every time step with ignition and low/high intenstiy burnt
  cat("... writing output layers", "\n")
  nfire <- sum(track.fire$year==t, na.rm=T)
  sizes <- filter(track.fire, year==t) %>% group_by(swc, fire.id) %>% summarise(ab=aburnt.highintens+aburnt.lowintens)
  # Ignitions' cell.id 
  igni.id <- burnt.cells[c(1,cumsum(sizes$ab)[1:(nfire-1)]+1)] 
  if(write.sp.outputs){
    BURNT[!is.na(MASK[])] <- land$distype*(land$tsdist==1)
    BURNT[igni.id] <- 9
    writeRaster(BURNT, paste0(out.path, "/lyr/DistType_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
  }
  
  save(land, file="./mdl_coupling/output/land.rdata")
  ## Deallocate memory
  gc()  
  cat("\n")
  
} 