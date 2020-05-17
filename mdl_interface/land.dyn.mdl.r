######################################################################################
##
######################################################################################

land.dyn.mdl <- function(scn.name){
  
  ## Load required packages and functions 
  suppressPackageStartupMessages({
    library(tictoc)
    library(sp)
    library(raster)  
    library(RANN)  # for nn2()
    library(Rcpp)
    library(tidyverse)
  })
  source("mdl/update.clim.r")
  source("mdl/update.interface.r")
  source("mdl/land.cover.change.r")
  source("mdl/prob.igni.r")
  source("mdl/growth.r")
  source("mdl/drought.r")
  source("mdl/cohort.establish.r")
  source("mdl/afforestation.r")
  source("mdl/forest.mgmt.r")
  source("mdl/fire.regime.r")
  source("mdl/post.fire.r")
  source("mdl/auxiliars.r")
  sourceCpp("mdl/is.in.cpp")
  
  ##To avoid library clashes
  select <- dplyr::select
  
  ## Load scenario definition (global variables and scenario parameters)
  ## and customized scenario parameters
  source(paste0("outputs/", scn.name, "/scn.def.r"))
  if(file.exists(paste0("outputs/", scn.name, "/scn.custom.def.r")))
    source(paste0("outputs/", scn.name, "/scn.custom.def.r"))
  
  
  ## Load:
  ## 1. Mask of the study area (raster)
  ## 2. Data frame with cell.id and coordinates x, y
  ## 3. Data frame of the model static variables 
  ## 4. Data frame with interface value
  load("inputlyrs/rdata/mask.rdata")
  load("inputlyrs/rdata/coordinates.rdata")
  load("inputlyrs/rdata/orography.rdata")
  # load("inputlyrs/rdata/harvest.rdata")
  load("inputlyrs/rdata/interface.rdata")
  
  
  ## Set the directory for writing spatial outputs (create it, if it does not exist yet) 
  if(write.sp.outputs){      
    if(!file.exists(paste0(out.path, "/lyr")))
      dir.create(file.path(getwd(), out.path, "/lyr"), showWarnings = F) 
  }

  
  ## List the name of the forest species
  species <- c("phalepensis", "pnigra", "ppinea", "psylvestris", "ppinaster", "puncinata",
               "aalba", "qilex", "qsuber", "qfaginea", "qhumilis", "fsylvatica", "other")
                  
  
  ## Translation equations from Basal Area to Volum, Volum with bark and Carbon
  eq.ba.vol <- read.table("inputfiles/EqBasalAreaVol.txt", header=T)
  eq.ba.volbark <- read.table("inputfiles/EqBasalAreaVolWithBark.txt", header=T)
  eq.ba.carbon <- read.table("inputfiles/EqBasalAreaCarbon.txt", header=T)
  
  
  ## Climatic severity and pctg hot days tabes
  clim.severity <- read.table(paste0("inputfiles/", file.clim.severity, ".txt"), header=T)
  
  
  ## Build the baseline time sequence and the time sequence of the processes (shared for all runs). 
  ## 1. Climate change, 2. Land-cover changes, 3. Forest management
  ## 4. Wildfires, 5. Prescribed burns, 6. Drought, 7. Post-fire regeneration,
  ## 8. Cohort establihsment, 9. Afforestation, 10. Growth
  time.seq <- seq(1, time.horizon, 1)
  if(time.horizon==1)
    clim.schedule <- 1
  else
    clim.schedule <- seq(1, time.horizon-1, clim.step)
  lchg.schedule <- seq(1, time.horizon, lchg.step)
  fmgmt.schedule <- seq(1, time.horizon, fmgmt.step)
  fire.schedule <- seq(1, time.horizon, fire.step)
  pb.schedule <- seq(1, time.horizon, pb.step)
  drought.schedule <- seq(1, time.horizon, drought.step)
  post.fire.schedule <- seq(1, time.horizon, post.fire.step)
  cohort.schedule <- seq(1, time.horizon, cohort.step)
  afforest.schedule <- seq(1, time.horizon, afforest.step)
  growth.schedule <- seq(1, time.horizon, growth.step)
  

  ## Tracking data.frames
  track.fmgmt <- data.frame(run=NA, year=NA, spp=NA, sylvi=NA, sawlog=NA, wood=NA)
  track.fire <-  data.frame(run=NA, year=NA, swc=NA, clim.sever=NA, fire.id=NA, fst=NA, wind=NA, atarget=NA, 
                            aburnt.highintens=NA, aburnt.lowintens=NA, asupp.fuel=NA, asupp.sprd=NA)
  track.fire.spp <-  data.frame(run=NA, year=NA, fire.id=NA, spp=NA, aburnt=NA, bburnt=NA)
  track.pb <-  data.frame(run=NA, year=NA, clim.sever=NA, fire.id=NA, 
                          wind=NA, atarget=NA, aburnt.lowintens=NA)
  track.drougth <- data.frame(run=NA, year=NA, spp=NA, ha=NA)
  track.cohort <- data.frame(run=NA, year=NA, spp.out=NA, Var2=NA, Freq=NA)
  track.post.fire <- data.frame(run=NA, year=NA, spp.out=NA, Var2=NA, Freq=NA)
  track.afforest <- data.frame(run=NA, year=NA, Var1=NA, Freq=NA)
  track.land <- data.frame(run=NA, year=NA, spp=NA, area=NA, vol=NA, volbark=NA, carbon=NA)
  
  
  ## Start the simulations   
  irun=1   # for testing
  for(irun in 1:nrun){
    
    ## Copy the schedulings in auxiliar vectors (only for those processes included in the current version)
    temp.clim.schedule <- clim.schedule
    temp.lchg.schedule <- lchg.schedule
    temp.fmgmt.schedule <- fmgmt.schedule
    temp.fire.schedule <- fire.schedule
    temp.pb.schedule <- pb.schedule
    temp.drought.schedule <- drought.schedule
    temp.post.fire.schedule <- post.fire.schedule
    temp.cohort.schedule <- cohort.schedule
    temp.afforest.schedule <- afforest.schedule
    temp.growth.schedule <- growth.schedule
    
    
    ## Load initial spatial dynamic state variables in a data.frame format
    load("inputlyrs/rdata/land.rdata")
    
    
    ## Start the discrete time sequence 
    for(t in time.seq){
      
      ## Track scenario, replicate and time step
      cat(paste0("scn: ", scn.name," - run: ", irun, "/", nrun, " - time: ", t, "/", time.horizon), "\n")
      
      
      ## 1. CLIMATE CHANGE  
      if(processes[clim.id] & t %in% temp.clim.schedule){
        clim <- update.clim(MASK, land, orography, decade=(1+floor(t/10))*10, clim.scn, clim.mdl)
        load(paste0("inputlyrs/rdata/sdm_base_", clim.scn, "_", clim.mdl, "_", (1+floor(t/10))*10, ".rdata"))
        temp.clim.schedule <- temp.clim.schedule[-1] 
      }
      
      
      ## 2. LAND-COVER CHANGE
      if(processes[lchg.id] & t %in% temp.lchg.schedule){
        # tic("  t")
        # Urbanization
        chg.cells <- land.cover.change(land, coord, interface, 1, t, numeric())
        land$spp[land$cell.id %in% chg.cells] <- 20 # urban
        land$biom[land$cell.id %in% chg.cells] <- NA
        land$age[land$cell.id %in% chg.cells] <- NA
        land$tsdist[land$cell.id %in% visit.cells] <- NA  # don't care the time since it's urban
        land$distype[land$cell.id %in% chg.cells] <- lchg.urb
        land$tburnt[land$cell.id %in% chg.cells] <- NA
        # Agriculture conversion
        visit.cells <- chg.cells
        chg.cells <- land.cover.change(land, coord, interface, 2, t, visit.cells)
        land$spp[land$cell.id %in% chg.cells] <- 16 # arableland or 17 - permanent crops
        land$biom[land$cell.id %in% chg.cells] <- NA
        land$age[land$cell.id %in% chg.cells] <- NA
        land$distype[land$cell.id %in% chg.cells] <- lchg.crp
        land$tsdist[land$cell.id %in% visit.cells] <- 0
        land$tburnt[land$cell.id %in% chg.cells] <- 0
        # Rural abandonment
        visit.cells <- c(visit.cells, chg.cells)
        chg.cells <- land.cover.change(land, coord, interface, 3, t, visit.cells)
        land$spp[land$cell.id %in% chg.cells] <- 14  # shrub
        land$biom[land$cell.id %in% chg.cells] <- 0
        land$age[land$cell.id %in% chg.cells] <- 0
        land$distype[land$cell.id %in% chg.cells] <- lchg.nat
        land$tsdist[land$cell.id %in% visit.cells] <- 0
        land$tburnt[land$cell.id %in% chg.cells] <- 0
        # Update interface values
        # toc()
        interface <- update.interface(land)
        temp.lchg.schedule <- temp.lchg.schedule[-1] 
      }
      
      
      ## 3. FOREST MANAGEMENT (under development)
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
      fintensity <- integer()
      fire.ids <- integer()
      id.fire <- annual.burnt <- 0
      if(processes[fire.id] & t %in% temp.fire.schedule){
        pigni <- prob.igni(land, orography, clim, interface)
        # Decide climatic severity of the year (default is mild)
        clim.sever <- 0
        if(runif(1,0,100) < clim.severity[clim.severity$year==t, ncol(clim.severity)]) # not-mild
          clim.sever <- 1
        # swc = wind, heat and regular. annual.burnt needed to compute PB target area 
        for(swc in 1:3){
          fire.out <- fire.regime(land, coord, orography, pigni, swc, clim.sever, t, 
                                  burnt.cells, fintensity, fire.ids, id.fire, annual.burnt)
          burnt.cells <- fire.out[[1]]; fintensity <- fire.out[[2]]; 
          fire.ids <- fire.out[[3]]; id.fire <- id.fire+nrow(fire.out[[4]])
          # track fire events and total annual burnt area
          if(nrow(fire.out[[4]])>0)
            track.fire <- rbind(track.fire, data.frame(run=irun, fire.out[[4]]))
          annual.burnt <- annual.burnt+sum(fire.out[[4]]$aburnt.highintens + fire.out[[4]]$aburnt.lowintens)
        }
        # track spp and biomass burnt
        aux <- data.frame(cell.id=burnt.cells, fire.id=fire.ids, fintensity) %>% 
               left_join(select(land, cell.id, spp, biom), by="cell.id") %>%
               mutate(bburnt=ifelse(fintensity>fire.intens.th, biom, biom*(1-fintensity))) %>%
               group_by(fire.id, spp) %>% summarize(aburnt=length(spp), bburnt=round(sum(bburnt, na.rm=T),1))
        if(nrow(aux)>0)
          track.fire.spp <-  rbind(track.fire.spp, data.frame(run=irun, year=t, aux)) 
        # Done with fires! When high-intensity fire, age = biom = 0 and dominant tree species may change
        # when low-intensity fire, age remains, spp remains and biomass.t = biomass.t-1 * (1-fintensity)
        burnt.intens <- fintensity>fire.intens.th
        land$tsdist[land$cell.id %in% burnt.cells] <- 0
        land$tburnt[land$cell.id %in% burnt.cells] <- land$tburnt[land$cell.id %in% burnt.cells] + 1
        land$distype[land$cell.id %in% burnt.cells[burnt.intens]] <- hfire
        land$distype[land$cell.id %in% burnt.cells[!burnt.intens]] <- lfire
        land$age[land$cell.id %in% burnt.cells[burnt.intens]] <- 0
        land$biom[land$cell.id %in% burnt.cells[burnt.intens]] <- 0
        land$biom[land$cell.id %in% burnt.cells[!burnt.intens]] <- 
          land$biom[land$cell.id %in% burnt.cells[!burnt.intens]]*(1-fintensity[!burnt.intens])
        temp.fire.schedule <- temp.fire.schedule[-1] 
        rm(fire.out)
      }
      
      
      ## 5. PRESCRIBED BURNS
      if(processes[pb.id] & t %in% temp.pb.schedule){
        fire.out <- fire.regime(land, coord, orography, pigni, swc=4, clim.sever, t, 
                                burnt.cells, burnt.intens, fintensity, fire.ids,  annual.burnt)
        pb.cells <- fire.out[[1]]
        if(nrow(fire.out[[4]])>0)
          track.pb <- rbind(track.pb, data.frame(run=irun, fire.out[[4]][,c(1,3,4,6,7,9)]))
        # done with prescribed burns!
        pb.cells <- pb.cells[!(pb.cells %in% burnt.cells)]
        land$tsdist[land$cell.id %in% pb.cells] <- 0
        land$tburnt[land$cell.id %in% pb.cells] <- land$tburnt[land$cell.id %in% pb.cells] + 1
        land$distype[land$cell.id %in% pb.cells] <- pb
        land$biom[land$cell.id %in% burnt.cells[!as.logical(burnt.intens)]] <- 
          land$biom[land$cell.id %in% burnt.cells[!as.logical(burnt.intens)]]*(1-fintensity[!as.logical(burnt.intens)])
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
        ## forest transition of tree species burnt in high intensity
        aux  <- post.fire(land, coord, orography, clim, sdm)
        if(nrow(aux)>0){
          spp.out <- land$spp[land$cell.id %in% aux$cell.id]
          land$spp[land$cell.id %in% aux$cell.id] <- aux$spp
          clim$spp[clim$cell.id %in% aux$cell.id] <- aux$spp
          clim$sdm[clim$cell.id %in% aux$cell.id] <- 1
          clim$sqi[clim$cell.id %in% aux$cell.id] <- aux$sqi
          track.post.fire <- rbind(track.post.fire, data.frame(run=irun, year=t, table(spp.out, aux$spp)))  
        }
        rm(aux) 
        temp.post.fire.schedule <- temp.post.fire.schedule[-1] 
      }
      
      
      ## 8. COHORT ESTABLISHMENT
      if(processes[cohort.id] & t %in% temp.cohort.schedule & length(killed.cells)>0){
        aux  <- cohort.establish(land, coord, orography, clim, sdm)
        spp.out <- land$spp[land$cell.id %in% killed.cells]
        land$spp[land$cell.id %in% killed.cells] <- aux$spp
        land$age[land$cell.id %in% aux$cell.id] <- 0  ## not sure if 0 or decade - t -1
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
        land$age[land$cell.id %in% aux$cell.id] <- 0
        land$tsdist[land$cell.id %in% aux$cell.id] <- 0
        land$distype[land$cell.id %in% aux$cell.id] <- afforest
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
      if(write.sp.outputs){
        MAP <- MASK
        cat("... writing output layers", "\n")
        nfire <- sum(track.fire$year==t, na.rm=T)
        sizes <- filter(track.fire, year==t) %>% group_by(swc, fire.id) %>% summarise(ab=aburnt.highintens+aburnt.lowintens)
        # Ignitions' cell.id 
        igni.id <- burnt.cells[c(1,cumsum(sizes$ab)[1:(nfire-1)]+1)] 
        MAP[!is.na(MASK[])] <- land$distype*(land$tsdist==1)
        MAP[igni.id] <- 9
        writeRaster(MAP, paste0(out.path, "/lyr/DistType_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
      }
      
      ## Deallocate memory
      gc()  
      cat("\n")
      
    } # time
  
    # Print maps at the end of the simulation period per each run
    if(write.sp.outputs){
      MAP <- MASK
      MAP[!is.na(MASK[])] <- land$tburnt
      writeRaster(MAP, paste0(out.path, "/lyr/TimesBurnt_r", irun, ".tif"), format="GTiff", overwrite=T)
    }

  } # run
  
  cat("... writing outputs", "\n")
  write.table(track.fmgmt[-1,], paste0(out.path, "/Management.txt"), quote=F, row.names=F, sep="\t")
  track.fire$rem <- pmax(0,track.fire$atarget-track.fire$aburnt.highintens-track.fire$aburnt.lowintens)
  write.table(track.fire[-1,], paste0(out.path, "/Fires.txt"), quote=F, row.names=F, sep="\t")
  write.table(track.fire.spp[-1,], paste0(out.path, "/FiresSpp.txt"), quote=F, row.names=F, sep="\t")
  write.table(track.pb[-1,], paste0(out.path, "/PrescribedBurns.txt"), quote=F, row.names=F, sep="\t")
  write.table(track.drougth[-1,], paste0(out.path, "/Drought.txt"), quote=F, row.names=F, sep="\t")
  names(track.post.fire)[4:5] <- c("spp.in", "ha")
  write.table(track.post.fire[-1,], paste0(out.path, "/PostFire.txt"), quote=F, row.names=F, sep="\t")
  names(track.cohort)[4:5] <- c("spp.in", "ha")
  write.table(track.cohort[-1,], paste0(out.path, "/Cohort.txt"), quote=F, row.names=F, sep="\t")
  names(track.afforest)[3:4] <- c("spp", "ha")
  write.table(track.afforest[-1,], paste0(out.path, "/Afforestation.txt"), quote=F, row.names=F, sep="\t")
  write.table(track.land[-1,], paste0(out.path, "/Land.txt"), quote=F, row.names=F, sep="\t")

}