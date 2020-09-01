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

  ##load global parameters (climate scenario and modules to be active)
  source("mdl_interface/global_parameters.R")

  if (MEDFIRE){
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
  }

  ##IPM functions
  if (IPM) {
	  source("./IPM/auxiliary_functions_v10 - Roberto.R")
	  source("./IPM/IPM_functions_v20_old.R")
	  source("./mdl_interface/build.var.IPM.R")
	  source("./mdl_interface/read_IPM_age.R")
	  sourceCpp("./IPM/IPM_functions_v25_1_Year_par.cpp")
	  ##IPM global variables and parameters
	  source("./mdl_interface/IPM_parameters.R")
  }
  
  ##To avoid library clashes
  select <- dplyr::select
	  
  time.seq <- seq(-9, time.horizon-10, 1) # From -9 to 0 runs IPM from 2000 to 2009
  
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

	  ##Tables to update sqi for dyn change of land spp
	  site.quality.spp <- read.table("Medfire/inputfiles/SiteQualitySpp.txt", header=T)
 	  site.quality.index <- read.table("Medfire/inputfiles/SiteQualityIndex.txt", header=T)
  	  site.quality.shrub <- read.table("Medfire/inputfiles/SiteQualityShrub.txt", header=T)
	  
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
  	num.trees.df <- ba
  	num.trees.df[,]  <- 0

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
		 
		 NUM_PLOTS <- 2
		 target<-c(1,2)
		 map <- map[target,]
		 ba <- ba[target,]
		 saplings <- saplings[target,]
		 IPM.forest.age <- IPM.forest.age[target,]
		 adult.trees<- lapply(adult.trees, function(x) {x[target,]})
	  }

	  if(change.IPM.plots.from.Medfire.LCT){
	    if (!MEDFIRE){ 
	      load("Medfire/inputlyrs/rdata/land.rdata")
	    }  
	  	empty.plots.IPM.index <- which(apply(ba,1,sum)==0)
	  	for(i in empty.plots.IPM.index){
	  		Medfire.spp <- land$spp[land$cell.id == map$Medfire.id[i]]
	  		if (Medfire.spp!=14){
	  			IPM.spp <- Medfire.index.IPM.spp[Medfire.spp]
	  			saplings[i,IPM.spp] <- new.saplings[IPM.spp]
	  			IPM.forest.age[i,IPM.spp]<- 10
	  		}
	  	}
	  	cat("empty plots filled according to MEDFIRE spp\n")

	  	tree.plots.to.delete.IPM.index <- which(map$Medfire.id %in% land$cell.id[land$spp==14 & land$cell.id %in% map$Medfire.id])
  		ba[tree.plots.to.delete.IPM.index,] <- 0
  		saplings[tree.plots.to.delete.IPM.index,] <- 0
  		for (s in 1:NUM_SP){
            adult.trees[[s]][tree.plots.to.delete.IPM.index,]<-0
          } #for species 
      IPM.forest.age[tree.plots.to.delete.IPM.index,] <- 0
      cat("tree plots that are shrub MEDFIRE LCT, destroyed spp\n")
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
    tic()
    ## Start the discrete time sequence 
    for(t in time.seq){
      
      ## Track scenario, replicate and time step
      cat(paste0("scn: ", scn.name," - run: ", irun, "/", nrun, " - time: ", t, "/", time.horizon), "\n")

      ## 0. Update land LCT with IPM dominant species at year 2010
      if(MEDFIRE & IPM.LCT & iyear == 2010) {
      	for( i in 1:nrow(map)){
      		if(sum(ba[i,])!=0){
      			land$spp[land$cell.id == map$Medfire.id[i]] <- IPM.index.Medfire.spp[which.max(ba[i,])]
      		}
      		else if(sum(saplings[i,])!=0){
      			land$spp[land$cell.id == map$Medfire.id[i]] <- IPM.index.Medfire.spp[which.max(saplings[i,])]
      		}
      		else{
      			land$spp[land$cell.id == map$Medfire.id[i]] <- 14
      		}
      	}
        cat("Land spp initialized to be same as IPM dominant spp\n")
      }
      
      
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
          if(iyear < 2019 & burn.hist.fires){
          	cat("Analysing historic burnts \n")
          	hist_fires <- raster(paste0("./Medfire/historic_fires/Fires_",iyear,".tif"))
          	burnt.cells <- which(!is.na(hist_fires[]))
          	burnt.intens <- rep(T, length(burnt.cells)) ##assumes all fires have high intensity but it should be cheked!
          	if (iyear>=2010){
          		##Only change Medfire land variable if year >=2010
	            land$age[land$cell.id %in% burnt.cells[burnt.intens]] <- 0
	            land$biom[land$cell.id %in% burnt.cells[burnt.intens]] <- 0
	            land$tsdist[land$cell.id %in% burnt.cells] <- 0
	          	land$tburnt[land$cell.id %in% burnt.cells] <- land$tburnt[land$cell.id %in% burnt.cells] + 1
	          	land$distype[land$cell.id %in% burnt.cells] <- hfire
	          	#land$distype[land$cell.id %in% burnt.cells[!burnt.intens]] <- lfire
        	}
          } else if (iyear>=2010){
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
	          land$biom[land$cell.id %in% burnt.cells[burnt.intens]] <- 0
	          land$biom[land$cell.id %in% burnt.cells[!burnt.intens]] <- 
	            land$biom[land$cell.id %in% burnt.cells[!burnt.intens]]*(1-fintensity[!burnt.intens])

              rm(fire.out); rm(aux)
          }
          temp.fire.schedule <- temp.fire.schedule[-1] 
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
                if (IPM.forest.age[i,j]>=reg.age[j] & fire.regeneration[j] ){  #if species can regenerate
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
            IPM.forest.age[burnt.cells.IPM.index,]<-0
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
            # Reset age of cells burnt in high intensity
        	land$age[land$cell.id %in% burnt.cells[burnt.intens] & !is.na(land$spp) & land$spp<=14] <- 0
            temp.post.fire.schedule <- temp.post.fire.schedule[-1] 
        	rm(aux); rm(spp.out)
        	if(IPM.post.fire){
        		## LCT is updated yearly for cell that have been burnt at least once
        		## Crop lands are excluded since they are not in IPM
        		IPM.burnt.plots.indexes <- which(map$Medfire.id %in% land$cell.id[land$tburnt>0 & !is.na(land$tburnt)])
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
        		##update sqi
        		Medfire.burnt.plots <- data.frame(cell.id=land$cell.id[land$cell.id %in%  map$Medfire.id[IPM.burnt.plots.indexes]],
        			         						spp=land$spp[land$cell.id %in%  map$Medfire.id[IPM.burnt.plots.indexes]])
        		Medfire.burnt.plots <- left_join(Medfire.burnt.plots, select(clim, cell.id, temp, precip), by = "cell.id") %>% 
						      left_join(select(orography, cell.id, aspect, slope), by = "cell.id") %>%
						      left_join(site.quality.spp, by = "spp") %>% left_join(site.quality.index, by = "spp") %>% 
						      mutate(aux=c0+c_mnan*temp+c2_mnan*temp*temp+c_plan*precip+c2_plan*precip*precip+c_aspect*ifelse(aspect!=1,0,1)+c_slope*slope/10) %>%
						      mutate(sq=1/(1+exp(-1*aux))) %>% mutate(sqi=ifelse(sq<=p50, 1, ifelse(sq<=p90, 2, 3))) %>%
						      select(cell.id, spp, temp, precip, sqi)
				sqi.shrub <- filter(Medfire.burnt.plots, spp==14) %>% select(spp, temp, precip) %>% left_join(site.quality.shrub, by = "spp") %>%
			      mutate(aux.brolla=c0_brolla+c_temp_brolla*temp+c_temp2_brolla*temp*temp+c_precip_brolla*precip+c_precip2_brolla*precip*precip,
			             aux.maquia=c0_maquia+c_temp_maquia*temp+c_temp2_maquia*temp*temp+c_precip_maquia*precip+c_precip2_maquia*precip*precip,
			             aux.boix=c0_boix+c_temp_boix*temp+c_temp2_boix*temp*temp+c_precip_boix*precip+c_precip2_boix*precip*precip,
			             sq.brolla=1/(1+exp(-1*aux.brolla)), sq.maquia=1/(1+exp(-1*aux.maquia)), sq.boix=1/(1+exp(-1*aux.boix)),
			             sqi=ifelse(sq.brolla>=sq.maquia & sq.brolla>=sq.maquia, 1,
			                        ifelse(sq.maquia>=sq.brolla & sq.maquia>=sq.boix, 2,
			                               ifelse(sq.boix>=sq.brolla & sq.boix>=sq.maquia, 3, 0))) )
			    Medfire.burnt.plots$sqi[Medfire.burnt.plots$spp==14] <- sqi.shrub$sqi
			    clim$sqi[clim$cell.id %in% map$Medfire.id[IPM.burnt.plots.indexes]] <- Medfire.burnt.plots$sqi
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
      if(IPM) {
        if(COLONIZATION & iyear>=2010){
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
                if(i %in% conifers) my.model <- colonization.glm[[1]]
                else if(i %in% quercus) my.model <- colonization.glm[[2]]
                else if(i %in% deciduous) my.model <- colonization.glm[[3]]
                suitable$colonized <- predict(my.model,newdata = suitable,type = "response")
                suitable$colonized <- ifelse(suitable$colonized > colonization.threshold,1,0)
              #  suitable<- suitable[suitable$colonized==1,]
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
          colonized.plots.ID <- unique(colonized.plots.ID) ##There can be duplicates since a plot can be colonized by more than 1 spp
          save(colonized.plots.ID, file=paste0("./mdl_interface/output/",scn.name,"/colonized.plots.ID_",iyear,".rdata"))
        }#if colonization module is active
      }#if IPM
      if(MEDFIRE){
        if(processes[afforest.id] & t %in% temp.afforest.schedule){
            aux  <- afforestation(land, coord, orography, clim, sdm)
            if(length(aux)!=0){
              if (IPM & COLONIZATION & IPM.afforestation){
              	aux <- aux[!(aux$cell.id %in% map$Medfire.id),] ##only colonize cells outside of IPM study area
              }
              if (nrow(aux)>0){
                land$spp[land$cell.id %in% aux$cell.id] <- aux$spp
                land$age[land$cell.id %in% aux$cell.id] <- 0
                land$tsdist[land$cell.id %in% aux$cell.id] <- 0
                land$distype[land$cell.id %in% aux$cell.id] <- afforest
                clim$spp[clim$cell.id %in% aux$cell.id] <- aux$spp
                clim$sdm[clim$cell.id %in% aux$cell.id] <- 1
                clim$sqi[clim$cell.id %in% aux$cell.id] <- aux$sqi
                track.afforest <- rbind(track.afforest, data.frame(run=irun, year=t, table(aux$spp)))
              }
            }  
            temp.afforest.schedule <- temp.afforest.schedule[-1] 
	        if(IPM.afforestation){
	        	##Update state values of plots colonized by IPM
	          land$age[(land$cell.id %in% map$Medfire.id[map$ID %in% colonized.plots.ID]) & land$spp==14] <- 0
	          land$tsdist[land$cell.id %in% map$Medfire.id[map$ID %in% colonized.plots.ID]] <- 0
	          land$distype[land$cell.id %in% map$Medfire.id[map$ID %in% colonized.plots.ID]] <- afforest
	          clim$sdm[land$cell.id %in% map$Medfire.id[map$ID %in% colonized.plots.ID] & land$spp==14] <- 1
      			## Track LCT evolution of all the cells that have IPM plots that have been colonized at some point
	      		IPM.colonized.plots.indexes <- which(map$Medfire.id %in% land$cell.id[land$distype==afforest & !is.na(land$distype)])
	        		for (col.plot in IPM.colonized.plots.indexes) {
	        			if (sum(ba[col.plot,])>0){
	        				IPM.spp <- which.max(ba[col.plot,])
	        				Medfire.spp <- IPM.index.Medfire.spp[IPM.spp]
	        				land$spp[land$cell.id %in%  map$Medfire.id[col.plot]] <- Medfire.spp
	        				clim$spp[clim$cell.id %in%  map$Medfire.id[col.plot]] <- Medfire.spp
	        			} else if(sum(saplings[col.plot,])>0){
	        				IPM.spp <- which.max(saplings[col.plot,])
	        				Medfire.spp <- IPM.index.Medfire.spp[IPM.spp]
	        				land$spp[land$cell.id %in%  map$Medfire.id[col.plot]] <- Medfire.spp
	        				clim$spp[clim$cell.id %in%  map$Medfire.id[col.plot]] <- Medfire.spp
	        			} else{
	        				land$spp[land$cell.id %in%  map$Medfire.id[col.plot]] <- 14
	        				clim$spp[clim$cell.id %in%  map$Medfire.id[col.plot]] <- 14
	        			}
	        		}#for colonized plot
	        	##Update sqi
	        	Medfire.col.plots <- data.frame(cell.id=land$cell.id[land$cell.id %in%  map$Medfire.id[IPM.colonized.plots.indexes]],
        			         						spp=land$spp[land$cell.id %in%  map$Medfire.id[IPM.colonized.plots.indexes]])
        		Medfire.col.plots <- left_join(Medfire.col.plots, select(clim, cell.id, temp, precip), by = "cell.id") %>% 
						      left_join(select(orography, cell.id, aspect, slope), by = "cell.id") %>%
						      left_join(site.quality.spp, by = "spp") %>% left_join(site.quality.index, by = "spp") %>% 
						      mutate(aux=c0+c_mnan*temp+c2_mnan*temp*temp+c_plan*precip+c2_plan*precip*precip+c_aspect*ifelse(aspect!=1,0,1)+c_slope*slope/10) %>%
						      mutate(sq=1/(1+exp(-1*aux))) %>% mutate(sqi=ifelse(sq<=p50, 1, ifelse(sq<=p90, 2, 3))) %>%
						      select(cell.id, spp, temp, precip, sqi)
				sqi.shrub <- filter(Medfire.col.plots, spp==14) %>% select(spp, temp, precip) %>% left_join(site.quality.shrub, by = "spp") %>%
			      mutate(aux.brolla=c0_brolla+c_temp_brolla*temp+c_temp2_brolla*temp*temp+c_precip_brolla*precip+c_precip2_brolla*precip*precip,
			             aux.maquia=c0_maquia+c_temp_maquia*temp+c_temp2_maquia*temp*temp+c_precip_maquia*precip+c_precip2_maquia*precip*precip,
			             aux.boix=c0_boix+c_temp_boix*temp+c_temp2_boix*temp*temp+c_precip_boix*precip+c_precip2_boix*precip*precip,
			             sq.brolla=1/(1+exp(-1*aux.brolla)), sq.maquia=1/(1+exp(-1*aux.maquia)), sq.boix=1/(1+exp(-1*aux.boix)),
			             sqi=ifelse(sq.brolla>=sq.maquia & sq.brolla>=sq.maquia, 1,
			                        ifelse(sq.maquia>=sq.brolla & sq.maquia>=sq.boix, 2,
			                               ifelse(sq.boix>=sq.brolla & sq.boix>=sq.maquia, 3, 0))) )
			    Medfire.col.plots$sqi[Medfire.col.plots$spp==14] <- sqi.shrub$sqi
			    clim$sqi[clim$cell.id %in% map$Medfire.id[IPM.colonized.plots.indexes]] <- Medfire.col.plots$sqi
	      	}
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
          rm(aux); rm(aux.shrub)
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
                  if (iyear>=2009){
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
                      #dummy <- dgamma(y[,j],1/10, rate=param.ingrowth1[j])*(param.ingrowth2[j]/10)
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
                num.trees.df[i,j] <-0
                if (saplings[i,j]>0){
                  IPM.forest.age[i,j] <- IPM.forest.age[i,j] + 1 
                }
              }else{
                ba[i,j] <- quadTrapezCpp_1(dummy*x2[,j],h[j],nx)*pi
                num.trees.df[i,j] <- num.trees
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
        cat(paste0("Saving IPM variables year:", iyear, "\n"))
        adult.trees.file <- paste0("./mdl_interface/output/",scn.name,"/trees_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
        ba.file <- paste0("./mdl_interface/output/",scn.name,"/ba_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
        saplings.file <- paste0("./mdl_interface/output/", scn.name, "/saplings_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
        age.file <- paste0("./mdl_interface/output/", scn.name, "/age_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
        num.trees.file <- paste0("./mdl_interface/output/", scn.name, "/num_trees_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
        if (iyear%%10==0){save(adult.trees, file=adult.trees.file)}
        #save(adult.trees, file=adult.trees.file)
        save(ba, file=ba.file)
        save(saplings, file=saplings.file)
        save(IPM.forest.age, file=age.file)
        save(num.trees.df, file=num.trees.file)
        }
      }
      if (MEDFIRE){
        land.file <- paste0("./mdl_interface/output/",scn.name,"/land_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
        save(land, file=land.file)
      } 
      # if(MEDFIRE){
      #   ## Print maps every time step with ignition and low/high intenstiy burnt
      #   if(write.sp.outputs){
      #     MAP <- MASK
      #     cat("... writing output layers", "\n")
      #     nfire <- sum(track.fire$year==t, na.rm=T)
      #     sizes <- filter(track.fire, year==t) %>% group_by(swc, fire.id) %>% summarise(ab=aburnt.highintens+aburnt.lowintens)
      #     # Ignitions' cell.id 
      #     igni.id <- burnt.cells[c(1,cumsum(sizes$ab)[1:(nfire-1)]+1)] 
      #     MAP[!is.na(MASK[])] <- land$distype*(land$tsdist==1)
      #     MAP[igni.id] <- 9
      #     writeRaster(MAP, paste0(out.path, "/lyr/DistType_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
      #   }
      # }

      iyear <- iyear +1
      ## Deallocate memory
      gc()  
      cat("\n")
      
    } # time
    toc()
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