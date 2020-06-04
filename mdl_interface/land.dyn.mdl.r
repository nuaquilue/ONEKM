######################################################################################
##
######################################################################################

land.dyn.mdl <- function(scn.name){
  
  ## Load required packages and functions 
  suppressPackageStartupMessages({
    library(raster)  
    library(RANN)  
    library(Rcpp)
    library(tidyverse)
    library(ff) ## IPM
  })
  options(fftempdir = "./")

  ##Medfire functions
  source("Medfire/mdl/update.clim.r")
  source("Medfire/mdl/update.interface.r")
  source("Medfire/mdl/land.cover.change.r")
  source("Medfire/mdl/prob.igni.r")
  source("Medfire/mdl/growth.r")
  source("Medfire/mdl/drought.r")
  source("Medfire/mdl/cohort.establish.r")
  source("Medfire/mdl/afforestation.r")
  source("Medfire/mdl/forest.mgmt.r")
  source("Medfire/mdl/fire.regime.r")
  source("Medfire/mdl/post.fire.r")
  source("Medfire/mdl/auxiliars.r")
  sourceCpp("Medfire/mdl/is.in.cpp")

  ##IPM functions
  source("./IPM/auxiliary_functions_v10 - Roberto.R")
  source("./IPM/IPM_functions_v20_old.R")
  source("./mdl_interface/build.var.IPM.r")
  source("./mdl_interface/read_IPM_age.r")
  sourceCpp("./IPM/IPM_functions_v23_1_Year.cpp")
  
  ##To avoid library clashes
  select <- dplyr::select
  
  ##load global parameters (climate scenario and modules to be active)
  source("./mdl_interface/global_parameters.r")

  if(IPM){
	  ##IPM global variables and parameters
	  source(paste0("mdl_interface/IPM_parameters.r"))
  }
	  
  time.seq <- seq(-9, time.horizon, 1) # From -9 to 0 runs IPM from 2000 to 2009
  
  if(MEDFIRE){
	  ## Load scenario definition (global variables and scenario parameters)
	  ## and customized scenario parameters
	  source(paste0("Medfire/outputs/", scn.name, "/scn.def.r"))
	  if(file.exists(paste0("Medfire/outputs/", scn.name, "/scn.custom.def.r")))
	    source(paste0("Medfire/outputs/", scn.name, "/scn.custom.def.r"))
	  ## Load Medfire variables:
	  ## 1. Mask of the study area (raster)
	  ## 2. Data frame with cell.id and coordinates x, y
	  ## 3. Data frame of the model static variables 
	  ## 4. Data frame with interface value
	  load("Medfire/inputlyrs/rdata/mask.rdata")
	  load("Medfire/inputlyrs/rdata/coordinates.rdata")
	  load("Medfire/inputlyrs/rdata/orography.rdata")
	  # load("Medfire/inputlyrs/rdata/harvest.rdata")
	  load("Medfire/inputlyrs/rdata/interface.rdata")
	  ## List the name of the forest species
  	  species <- c("phalepensis", "pnigra", "ppinea", "psylvestris", "ppinaster", "puncinata",
               "aalba", "qilex", "qsuber", "qfaginea", "qhumilis", "fsylvatica", "other")
  	  ## Translation equations from Basal Area to Volum, Volum with bark and Carbon
	  eq.ba.vol <- read.table("Medfire/inputfiles/EqBasalAreaVol.txt", header=T)
	  eq.ba.volbark <- read.table("Medfire/inputfiles/EqBasalAreaVolWithBark.txt", header=T)
	  eq.ba.carbon <- read.table("Medfire/inputfiles/EqBasalAreaCarbon.txt", header=T)
	  
	  ## Climatic severity and pctg hot days tabes
	  clim.severity <- read.table(paste0("Medfire/inputfiles/", file.clim.severity, ".txt"), header=T)
	  
	  ## Build the baseline time sequence and the time sequence of the processes (shared for all runs). 
	  ## 1. Climate change, 2. Land-cover changes, 3. Forest management
	  ## 4. Wildfires, 5. Prescribed burns, 6. Drought, 7. Post-fire regeneration,
	  ## 8. Cohort establihsment, 9. Afforestation, 10. Growth
	  if(time.horizon==1)
	    clim.schedule <- 1
	  else
	    clim.schedule <- seq(1, time.horizon-1, clim.step)
	  lchg.schedule <- seq(1, time.horizon, lchg.step)
	  fmgmt.schedule <- seq(1, time.horizon, fmgmt.step)
	  fire.schedule <- seq(-9, time.horizon, fire.step) #burns IPM
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
	  
	  ## Set the directory for writing spatial outputs (create it, if it does not exist yet) 
	  if(write.sp.outputs){      
	    if(!file.exists(paste0(out.path, "/lyr")))
	      dir.create(file.path(getwd(), out.path, "/lyr"), showWarnings = F) 
	  }
  }
  
  ## Start the simulations   
  for(irun in 1:nrun){
    if(MEDFIRE){
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
   		load("Medfire/inputlyrs/rdata/land.rdata")
	}
  	
  	if(IPM){
      ## Load IPM dynamic variables
	  if (file.exists(orig.adult.trees.file) & file.exists(orig.ba.file) & file.exists(orig.saplings.file)){
		  load(orig.adult.trees.file); load(orig.ba.file); load(orig.saplings.file) 
		  cat("IPM initial variables loaded\n")
	  } else{
	      cat("building IPM dynamic variables\n")
	  	  build.var.IPM(clim.scn, clim.mdl)
	  	  load(orig.adult.trees.file); load(orig.ba.file); load(orig.saplings.file) 
	      cat("IPM initial variables loaded\n")
	  }

	  ##IPM plots age
	  if (file.exists(orig.plots.age.file)){
	  	load(orig.plots.age.file) ## to do
	  } else{
	  	read.IPM.age(ba, map)
	    cat("building age variable\n")
	  	load(orig.plots.age.file)
	  }
    cat("IPM age loaded\n")
	  if(testing){
	  	 ##select a smaller number of cells within Barcelona
		 ini_cells <- c()
		 NUM_PLOTS <- 100
		 i <- 1
		 while(length(ini_cells)<NUM_PLOTS){
		   if (sum(ba[i,])>0){
		     ini_cells <- c(ini_cells, i)
		   }
		   i <- i +1
		 }
		 map <- map[ini_cells,]
		 ba <- ba[ini_cells,]
		 saplings <- saplings[ini_cells,]
		 IPM.forest.age <- IPM.forest.age[ini_cells,]
		 adult.trees<- lapply(adult.trees, function(x) {x[ini_cells,]})
		 
		 # NUM_PLOTS <- 11
		 # target<-c(9,10,25,30,77,79,81,90,93,94,96)
		 # map <- map[target,]
		 # ba <- ba[target,]
		 # saplings <- saplings[target,]
		 # IPM.forest.age <- IPM.forest.age[target,]
		 # adult.trees<- lapply(adult.trees, function(x) {x[target,]})
	  }

    
	  if(ALL.RESULTS){
	    results <- list()
	    for(i in 1:NUM_SP){
	        results[[i]] <- data.frame(ID = integer(nrow(map)),
	                                   adults_2000 = integer(nrow(map)),
	                                   deaths_2000 = integer(nrow(map)),
	                                   new_adults_2000 = integer(nrow(map)),
	                                   saplings_2000 = integer(nrow(map)),
	                                   basal_area_2000 = numeric(nrow(map)),
	                                   #
	                                   adults_2010 = integer(nrow(map)),
	                                   deaths_2010 = integer(nrow(map)),
	                                   new_adults_2010 = integer(nrow(map)),
	                                   saplings_2010 = integer(nrow(map)),
	                                   basal_area_2010 = numeric(nrow(map)),
	                                   #
	                                   adults_2020 = integer(nrow(map)),
	                                   deaths_2020 = integer(nrow(map)),
	                                   new_adults_2020 = integer(nrow(map)),
	                                   saplings_2020 = integer(nrow(map)),
	                                   basal_area_2020 = numeric(nrow(map)),
	                                   #
	                                   adults_2030 = integer(nrow(map)),
	                                   deaths_2030 = integer(nrow(map)),
	                                   new_adults_2030 = integer(nrow(map)),
	                                   saplings_2030 = integer(nrow(map)),
	                                   basal_area_2030 = numeric(nrow(map)),
	                                   #
	                                   adults_2040 = integer(nrow(map)),
	                                   deaths_2040 = integer(nrow(map)),
	                                   new_adults_2040 = integer(nrow(map)),
	                                   saplings_2040 = integer(nrow(map)),
	                                   basal_area_2040 = numeric(nrow(map)),
	                                   #
	                                   adults_2050 = integer(nrow(map)),
	                                   deaths_2050 = integer(nrow(map)),
	                                   new_adults_2050 = integer(nrow(map)),
	                                   saplings_2050 = integer(nrow(map)),
	                                   basal_area_2050 = numeric(nrow(map)),
	                                   #
	                                   adults_2060 = integer(nrow(map)),
	                                   deaths_2060 = integer(nrow(map)),
	                                   new_adults_2060 = integer(nrow(map)),
	                                   saplings_2060 = integer(nrow(map)),
	                                   basal_area_2060 = numeric(nrow(map)),
	                                   #
	                                   adults_2070 = integer(nrow(map)),
	                                   deaths_2070 = integer(nrow(map)),
	                                   new_adults_2070 = integer(nrow(map)),
	                                   saplings_2070 = integer(nrow(map)),
	                                   basal_area_2070 = numeric(nrow(map)),
	                                   #
	                                   adults_2080 = integer(nrow(map)),
	                                   deaths_2080 = integer(nrow(map)),
	                                   new_adults_2080 = integer(nrow(map)),
	                                   saplings_2080 = integer(nrow(map)),
	                                   basal_area_2080 = numeric(nrow(map)),
	                                   #  
	                                   adults_2090 = integer(nrow(map)),
	                                   deaths_2090 = integer(nrow(map)),
	                                   new_adults_2090 = integer(nrow(map)),
	                                   saplings_2090 = integer(nrow(map)),
	                                   basal_area_2090 = numeric(nrow(map)),
	                                   #
	                                   temperature = numeric(nrow(map)),
	                                   precipitation = numeric(nrow(map)))
	        results[[i]]$ID <- map$ID
	    }
	  }# if ALL.RESULTS
    }#if IPM    
    
    
    iyear <- 2000
    ## Start the discrete time sequence 
    for(t in time.seq){
      
      ## Track scenario, replicate and time step
      cat(paste0("scn: ", scn.name," - run: ", irun, "/", nrun, " - time: ", t, "/", time.horizon), "\n")
      
      
      ## 1. CLIMATE CHANGE  
      if(MEDFIRE){
        if(processes[clim.id] & t %in% temp.clim.schedule){
          clim <- update.clim(MASK, land, orography, decade=(1+floor(t/10))*10, clim.scn, clim.mdl)
          load(paste0("Medfire/inputlyrs/rdata/sdm_base_", clim.scn, "_", clim.mdl, "_", (1+floor(t/10))*10, ".rdata"))
          temp.clim.schedule <- temp.clim.schedule[-1] 
        }
      }
      if(IPM){
        ##Load IPM climate (clima_IPM)
        decade <-  10*(floor((iyear-2000)/10))
        if(decade==0)
          load(paste0("IPM/clima/climate_", clim.scn, "_", clim.mdl, "_00.rdata"))
        else
          load(paste0("IPM/clima/climate_", clim.scn, "_", clim.mdl, "_", decade, ".rdata"))

        my.id <- match(map$ID,clima_IPM$ID)
      }


      ## 2. LAND-COVER CHANGE
      if(MEDFIRE){
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
      }
      
      
      ## 3. FOREST MANAGEMENT (under development)
      if(MEDFIRE){
        if(processes[fmgmt.id] & t %in% temp.fmgmt.schedule){
          aux <- forest.mgmt(land, coord, clim, harvest, t)
          land$tsdist[land$cell.id %in% aux$cell.id] <- 0
          land$distype[land$cell.id %in% aux$cell.id] <- fmgmt.id*10+aux$sylvi
          track.fmgmt <- rbind(track.fmgmt, 
                               data.frame(run=irun, year=t, 
                                          group_by(aux, spp, sylvi) %>% summarize(sawlog=sum(vol.sawlog), wood=sum(vol.wood))))
          temp.fmgmt.schedule <- temp.fmgmt.schedule[-1] 
        }
      }
      
      
      ## 4. FIRE
      if(MEDFIRE){
        ## Tracking variables to be re-initialized each time step
        ## Out of the "if(fires)" in case only prescribed burns are applied
        burnt.cells <- integer()
        fintensity <- integer()
        fire.ids <- integer()
        id.fire <- annual.burnt <- 0
        if(processes[fire.id] & t %in% temp.fire.schedule){
          if(iyear < 2019){
          	hist_fires <- raster(paste0("./Medfire/historic_fires/Fires_",iyear,".TIF"))
          	burnt.cells <- which(!is.na(hist_fires[]))
          	land$tsdist[land$cell.id %in% burnt.cells] <- 0
            land$tburnt[land$cell.id %in% burnt.cells] <- land$tburnt[land$cell.id %in% burnt.cells] + 1
            land$distype[land$cell.id %in% burnt.cells] <- hfire
            #land$distype[land$cell.id %in% burnt.cells[!burnt.intens]] <- lfire
            land$age[land$cell.id %in% burnt.cells[burnt.intens]] <- 0
            land$biom[land$cell.id %in% burnt.cells[burnt.intens]] <- 0
          } else{
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
          }
          temp.fire.schedule <- temp.fire.schedule[-1] 
          rm(fire.out)
        }
      }
      

      ## 5. PRESCRIBED BURNS
      if(MEDFIRE){
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
      }
      
      
      ## 6. DROUGHT
      if(MEDFIRE){
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
      }
      
      
      ## 7. POST-FIRE REGENERATION
      if(MEDFIRE){
        if(IPM){#IPM only happens if Medfire
          ##IPM plots burnt by medfire
          burnt.cells.IPM.index <- which(map$Medfire.id %in% burnt.cells[burnt.intens]) 
          if(length(burnt.cells.IPM.index)!=0){
            #load(adult.trees.file); load(ba.file); load(saplings.file); load(future.saplings.file); load(IPM.forest.age.file)
            ##Calculate future saplings
            for (i in burnt.cells.IPM.index){
              tot.saplings <- sum(saplings[i,])
              for(j in 1:NUM_SP){
                if (IPM.forest.age[i,j]>9 & fire.regeneration[j] ){  #if species can regenerate
                  if (sum(ba[i,])>0){
                    ## saplings: fixed number of new.saplings times the abundance proportion of the species in the plot (in ba) 
                    saplings[i,j] <- new.saplings[j]*(ba[i,j]/sum(ba[i,]))
                  }
                  else{
                    if(sum(saplings[i,])==0){ ##In principle it cannot happen, used to debug
                      stop("error: burnt cell with age older 9 without saplings or ba")
                    }
                    saplings[i,j] <- new.saplings[j]*(saplings[i,j]/tot.saplings)
                  }}
                else {saplings[i,j] <- 0} #species cannot regenerate
              } #for tree species
            } #for burnt IPM plot
            ##Burns IPM plots by assigning adult.trees, ba and saplings to 0
            for (s in 1:NUM_SP){
                adult.trees[[s]][burnt.cells.IPM.index,]<-0
              } #for species 
            ba[burnt.cells.IPM.index,]<-0
            IPM.forest.age[burnt.cells.IPM.index]<-0
          } # if  there are high intensity burnt cells
          ##low intensity burnt cells
          low.intensity.burnt <- burnt.cells[!burnt.intens]
          burnt.cells.IPM.index <- which(map$Medfire.id %in% low.intensity.burnt) 
          if(length(burnt.cells.IPM.index)!=0){
          	burnt.IPM.i <- match(map$Medfire.id[burnt.cells.IPM.index], low.intensity.burnt)
          	for (s in 1:NUM_SP){
                adult.trees[[s]][burnt.cells.IPM.index,]<- adult.trees[[s]][burnt.cells.IPM.index,]*(1-fintensity[!burnt.intens][burnt.IPM.i])
              } #for species 
            ba[burnt.cells.IPM.index,]<- ba[burnt.cells.IPM.index,]*(1-fintensity[!burnt.intens][burnt.IPM.i])
            saplings[burnt.cells.IPM.index,]<- saplings[burnt.cells.IPM.index,]*(1-fintensity[!burnt.intens][burnt.IPM.i])
			}
		}##if IPM
        ##MEDFIRE
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
        	if(IPM.post.fire){
        		## LCT is updated yearly for cell that have been burnt at least once
        		IPM.burnt.plots.indexes <- which(map$Medfire.id %in% land[land$tburnt>0,cell.id])
        		for (burnt.plot in IPM.burnt.plots.indexes) {
        			if (sum(ba[burnt.plot,])>0){
        				IPM.spp <- which.max(ba[burnt.plot,])
        				Medfire.spp <- IPM.index.Medfire.spp[IPM.spp]
        				land$spp[land$cell.id %in%  map$Medfire.id[burnt.plot]] <- Medfire.spp
        				clim$spp[clim$cell.id %in%  map$Medfire.id[burnt.plot]] <- Medfire.spp
        			} else if(sum(saplings[burnt.plot,])>0){
        				IPM.spp <- which.max(saplings[burnt.plot,])
        				Medfire.spp <- IPM.index.Medfire.spp[IPM.spp]
        				land$spp[land$cell.id %in%  map$Medfire.id[burnt.plot]] <- Medfire.spp
        				clim$spp[clim$cell.id %in%  map$Medfire.id[burnt.plot]] <- Medfire.spp
        			} else{
        				land$spp[land$cell.id %in%  map$Medfire.id[burnt.plot]] <- 14 #shrub
        				clim$spp[clim$cell.id %in%  map$Medfire.id[burnt.plot]] <- 14 #shrub
        			}
        		}#for burnt plot
        	}#if IPM.post.fire
        }#if process
      }#if Medfire
      

      ## 8. COHORT ESTABLISHMENT
      if(MEDFIRE){
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
      }
      

      ## 9. AFFORESTATION
      if(MEDFIRE){
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
      }

      if(IPM) {
        if(COLONIZATION){
          # keep a dataframe with the basal area for potential colonization
          temp.suitable <- data.frame(ID = map$ID,
                                      saplings = integer(nrow(map)),
                                      plot_basal_area = integer(nrow(map)),
                                      #neigh = integer(nrow(map)),
                                      sp_age = integer(nrow(map)),
                                      neigh_ba = integer(nrow(map)),
                                      suitable = logical(nrow(map)))
          temp.suitable$plot_basal_area <- rowSums(ba) #sapply(X=trees,FUN=function(x) sum(unlist(x$ba)))
          # default 
          temp.suitable$suitable <- T
          colonized.plots.ID <- c()
          order <- sample(x = 1:NUM_SP,size = NUM_SP,replace = F)
          for(i in order){            
            print(paste(date()," - ",scn.name," - ",clim.scn," - year ",iyear," - colonization for sp ",i,"...",sep=""))
            temp.suitable$sp_age <- IPM.forest.age[,i]
            suitable <- GetSuitablePlots(map=map, temp.suitable=temp.suitable, ba = ba, distancias=distancias, sp=tesauro[i,1],
                                          max.dist = max.dist, verbose=T)   
            ##apply colonization
            if(!is.null(nrow(suitable))){
              if(nrow(suitable)>0){
              	if (MEDFIRE){
              		##Only suitable if it has been colonized by Medfire as well
              		if (processes[afforest.id]){
              			suitable <- left_join(suitable,map[,c("ID","Medfire.id")], by="ID") %>% 
              						left_join(land[,c("cell.id","spp")], by=c("Medfire.id"="cell.id"))
              			suitable <- suitable[(suitable$spp != 14),] #LCT that remain shrub after medfire colonization cannot be colonized
              		}
              	} ##if Medfire
                if(i %in% conifers) my.model <- colonization.glm[[1]]
                else if(i %in% quercus) my.model <- colonization.glm[[2]]
                else if(i %in% deciduous) my.model <- colonization.glm[[3]]
                suitable$colonized <- predict(my.model,newdata = suitable,type = "response")
                suitable$colonized <- ifelse(suitable$colonized > colonization.threshold,1,0)
                if (nrow(suitable)>0){
	                k <- match(suitable$ID,map$ID)
	                colonized.plots.ID <- c(colonized.plots.ID, map$ID[k])
	                #add new saplings to saplings list
	                for(j in 1:nrow(suitable)){
	                  saplings[k[j],tesauro[i,1]] <- new.saplings[tesauro[i,1]] 
	                }
            	}
              }
            }# if there are suitable plots
          }#for each species
          ##colonize plots colonized in Medfire but not in IPM
          colonized.plots.ID <- unique(colonized.plots.ID)
          if (MEDFIRE){
            if (processes[afforest.id]){
            	#colonized.medfire <- which(map$Medfire.id %in% aux$cell.id)
      			#to.colonize.ID <- map$ID[colonized.medfire][!(map$ID[colonized.medfire] %in% colonized.plots.ID )] 
      			to.colonize.IPM.i <- which(map$Medfire.id %in% aux$cell.id & !(map$ID %in% colonized.plots.ID )) 
      			if(length(to.colonize.IPM.i)>0){
      				#k <- match(to.colonize.ID,map$ID)
      				for(j in to.colonize.IPM.i){
      					spp.Medfire <- land$spp[land$cell.id==map$Medfire.id[to.colonize.IPM.i[j]]]
      					spp.IPM <- Medfire.index.IPM.spp[spp.Medfire]
      					saplings[to.colonize.IPM.i[j], spp.IPM] <- new.saplings[spp.IPM]
      				}
      			}
      		}
          } ## if Medfire
        }#if colonization module is active
      }#if IPM

      ## Track colonized plots LCT evolution of all the cells that have IPM plots
      if(MEDFIRE){
      	if(IPM.afforestation){
      		IPM.colonized.plots.indexes <- which(map$Medfire.id %in% land[land$distype==afforest,cell.id])
        		for (col.plot in IPM.colonized.plots.indexes) {
        			if (sum(ba[col.plot,])>0){
        				IPM.spp <- which.max(ba[burnt.plot,])
        				Medfire.spp <- IPM.index.Medfire.spp[IPM.spp]
        				land$spp[land$cell.id %in%  map$Medfire.id[col.plot]] <- Medfire.spp
        				clim$spp[clim$cell.id %in%  map$Medfire.id[col.plot]] <- Medfire.spp
        			} else if(sum(saplings[col.plot,])>0){
        				IPM.spp <- which.max(saplings[burnt.plot,])
        				Medfire.spp <- IPM.index.Medfire.spp[IPM.spp]
        				land$spp[land$cell.id %in%  map$Medfire.id[col.plot]] <- Medfire.spp
        				clim$spp[clim$cell.id %in%  map$Medfire.id[col.plot]] <- Medfire.spp
        			} else{
        				land$spp[land$cell.id %in%  map$Medfire.id[col.plot]] <- 14
        				clim$spp[clim$cell.id %in%  map$Medfire.id[col.plot]] <- 14
        			}
        		}#for colonized plot
      	}
      }          	
      
      ## 10. GROWTH
      if(MEDFIRE){
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
      }

      if (IPM){
        ## IPM: Growth, death, ingrowth and newsaplings
        if(ALL.RESULTS){
          # indexes for storing results
          my.new.adults <- which(names(results[[1]]) == paste("new_adults_",iyear,sep=""))
          my.adults <- which(names(results[[1]]) == paste("adults_",iyear,sep=""))
          my.deaths <- which(names(results[[1]]) == paste("deaths_",iyear,sep=""))
          my.saplings <- which(names(results[[1]]) == paste("saplings_",iyear,sep=""))
          my.ba <- which(names(results[[1]]) == paste("basal_area_",iyear,sep=""))
          previous.adults <- which(names(results[[1]]) == paste("adults_",(iyear - 10),sep=""))
        }

        ## IPM loop        
        for (i in 1:NUM_PLOTS) {
          
          if (round(i/10)*10==i) print(paste(date(),"- IPM: running loop",i,"of",NUM_PLOTS," - year",iyear,sep=" "))
          
          sapl <- saplings[i,] ##store old saplings so it can be used by ingrowth function
          
          ## are there adult trees or saplings?
          if (sum(ba[i,])>0 | sum(sapl)>0){
            
            temper <- clima_IPM$temp[my.id[i]]
            rain <- clima_IPM$precip[my.id[i]]
            anom.temper <- clima_IPM$anom_temp[my.id[i]]  
            anom.rain <- clima_IPM$anom_pl[my.id[i]]
            
            q <- c(rain,temper,sum(ba[i,]),anom.rain,anom.temper)

            ## if adult trees, apply Growth, Death and Saplings prediction
            if(sum(ba[i,])>0){
              
              sum.ntrees <- numeric(length(newspecies))
              
              for(i.sp in 1:length(newspecies)){
                ## this is not the number of trees, it's just for knowing if there are trees of that sp
                sum.ntrees[i.sp] <- sum(adult.trees[[i.sp]][i,])
              }
              
              param.survival2 <- CoefTimesVarCpp(param=survival.coef,z=q,type="survival")
              param.growth2 <- CoefTimesVarCpp(param=growth.coef,z=q,type="growth")
              param.sapl2 <- CoefTimesVarCpp(param=saplings.coef,z=q,spba=ba[i,],type="recruitment")
              
              # IPM, continuous case.
              order <- sample(x = 1:NUM_SP,size = NUM_SP,replace = F)
              for (j in order) {
                if (sum.ntrees[j]>0) {    # Are there trees of jth species?
                  dummy <- IPMadultSurvivalTimesGrowthCpp(ntrees=adult.trees[[j]][i,],#trees[[i]]$num[[j]],
                                                           x=x[,j],
                                                           y=y[,j],
                                                           param_survival=c(param.survival1[j],param.survival2[j]),
                                                           param_growth=c(param.growth1[j],param.growth2[j],param.growth3[j],param.growth4[j]),
                                                           t_diff=t.diff,
                                                           max_diam=max.diam[j],
                                                           h=h[j],
                                                           nx=nx)          
                  adult.trees[[j]][i,] <- dummy
                  
                  ##Calculate newsaplings 
                  if (iyear>2009){
                    if (IPM.forest.age[i,j]>29){
                      p.sapl <- c(param.sapl1[j],param.sapl2[j])
                      #careful with the saplings prediction. If not controlled, it can both explode and drop below 0
                      saplings[i,j] <- ifelse(IPMSaplingsCpp(sapl[j],p.sapl)>100,100,IPMSaplingsCpp(sapl[j],p.sapl))
                      if(saplings[i,j] < 0){
                        saplings[i,j] <- 0
                      }
                    }
                  }

                  if(ALL.RESULTS){
                    results[[j]][i,my.adults] <- quadTrapezCpp_1(dummy,h[j],nx)
                    results[[j]][i,my.saplings] <- ifelse(is.null(saplings[i,j]),
                                                          0,
                                                          saplings[i,j]*(10000/(pi*25)))
                    if(iyear == 2000){
                      results[[j]][i,my.deaths] <- ifelse(sum(trees.orig$factor[trees.orig$num.sp == j & trees.orig$ID == map$ID[i]]) - results[[j]]$adults_2000[i] > 0,
                                                          sum(trees.orig$factor[trees.orig$num.sp == j & trees.orig$ID == map$ID[i]]) - results[[j]]$adults_2000[i],
                                                          0)
                    }else{
                      results[[j]][i,my.deaths] <- ifelse(results[[j]][i,previous.adults] - results[[j]][i,my.adults] > 0,
                                                          results[[j]][i,previous.adults] - results[[j]][i,my.adults],
                                                          0)
                    }
                  }
                }#sum(ntrees[[j]])
              }# for j in NUM_SP
            }# if(sum(unlist(trees[[i]]$ba))>0)

            ## if saplings on the previous step, apply INGROWTH
            if(ingrowth){
              if(sum(sapl)>0){
  
                  param.ingrowth2 <-  CoefTimesVarCpp(param=ingrowth.coef,z=q,s=sapl,type="ingrowth")
  
                  order <- sample(x = 1:NUM_SP,size = NUM_SP,replace = F)
                  for(j in order){
  
                    if(sapl[j] > 0 & (!BASAL.AREA.THRESHOLD | (sum(ba[i,]) < BA_threshold$perc_95[j])) & IPM.forest.age[i,j]>9){
  
                      dummy <- IPMIngrowthIdentityCpp(y[,j],c(param.ingrowth1[j],param.ingrowth2[j]))
                      adult.trees[[j]][i,] <- adult.trees[[j]][i,] + dummy
  
                      if(ALL.RESULTS){
                        results[[j]][i,my.new.adults] <- quadTrapezCpp_1(dummy,h[j],nx)
                        results[[j]][i,my.adults] <- quadTrapezCpp_1(adult.trees[[j]][i,],h[j],nx)
                      }
                    }# if basal area threshold is met for j-th species and there are saplings older than 9 yo
                  } # for all j-th species
                } # if there are saplings of any sp
            }
            
            ## if any of the two conditions is met (saplings or adults), update basal area of the plot and age
            order <- sample(x = 1:NUM_SP,size = NUM_SP,replace = F)
            for (j in order) {
              dummy <- adult.trees[[j]][i,]
              num.trees <- quadTrapezCpp_1(dummy,h[j],nx)
              if(num.trees < 0.1){
                adult.trees[[j]][i,] <- 0
                ba[i,j] <- 0
                if (saplings[i,j]>0){
                  IPM.forest.age[i,j] <- IPM.forest.age[i,j] + 1 
                }
              }else{
                ba[i,j] <- quadTrapezCpp_1(dummy*x2[,j],h[j],nx)*pi
                IPM.forest.age[i,j] <- IPM.forest.age[i,j] + 1 
              }

              if(ALL.RESULTS){
                results[[j]][i,my.ba] <- ba[i,j] #ifelse(is.null(trees[[i]]$ba[[j]]),0,trees[[i]]$ba[[j]])
              }
            }# for j in NUM_SP
          }# if saplings or adults
          
        }# for i in plots
      }

      ##SAVE RESULTS
      if(IPM){
        if(save.IPM.variables){
        adult.trees.file <- paste0("./mdl_interface/output/trees_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
        ba.file <- paste0("./mdl_interface/output/ba_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
        saplings.file <- paste0("./mdl_interface/output/saplings_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
        save(adult.trees, file=adult.trees.file)
        save(ba, file=ba.file)
        save(saplings, file=saplings.file)
        #save(IPM.forest.age, file=IPM.forest.age.file)
        }
      }
      
      if(MEDFIRE){
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
      }

      iyear <- iyear +1
      ## Deallocate memory
      gc()  
      cat("\n")
      
    } # time
    
    if(MEDFIRE){
      # Print maps at the end of the simulation period per each run
      if(write.sp.outputs){
        MAP <- MASK
        MAP[!is.na(MASK[])] <- land$tburnt
        writeRaster(MAP, paste0(out.path, "/lyr/TimesBurnt_r", irun, ".tif"), format="GTiff", overwrite=T)
      }
    }

  } # run
  
  if(MEDFIRE){
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
}