######################################################################################
###  define.scenario()
###
###  Description >  Initialize the scenario parameters and the global variables of 
###   MEDFIRE model
###
###  Arguments >  
###   scn.name : identificative name of the scenario (string)
###
###  Details > By default, the output directory is ..\outputs\scn.name\ and all the
###   objects are saved in the file scn.def.r.
###
###  Value >  An R script.
######################################################################################

define.scenario <- function(scn.name){

  cat("Initializing parameters", "\n")

  ## Output directory (do not never change that, please!)
  out.path <- paste0("outputs/", scn.name)
  
  ## Flags to write spatial and tabular output data
  write.sp.outputs <- TRUE
  write.tbl.outputs <- TRUE
  
  ## Id for distrubance types
  lchg <- 1
  cut <- 2
  thin <- 3
  hfire <- 4
  lfire <- 5
  pb <- 6
  drght <- 7
  afforest <- 8
    
  ## Processes of the model included (TRUE-IN our FALSE-OUT),
  ## process ID, and frequency (in years) of each process
  processes <- c(TRUE,   # 1. Climate change
                 FALSE,  # 2. Land-cover changes
                 FALSE,  # 3. Forest management
                 TRUE,   # 4. Wildfires
                 FALSE,   # 5. Prescribed burns
                 FALSE,  # 6. Drought
                 FALSE,   # 7. Post-fire regeneration
                 FALSE,  # 8. Cohort establihsment
                 FALSE,  # 9. Afforestation
                 TRUE)   # 10. Growth
  clim.id <- 1; clim.step <- 10
  lchg.id <- 2; lchg.step <- 5
  fmgmt.id <- 3; fmgmt.step <- 1
  fire.id <- 4; fire.step <- 1
  pb.id <- 5; pb.step <- 1
  drought.id <- 6; drought.step <- 1
  post.fire.id <- 7; post.fire.step <- 1
  cohort.id <- 8; cohort.step <- 1
  afforest.id <- 9; afforest.step <- 1
  growth.id <- 10; growth.step <- 1
  
  ## Time lenght (in years) of a model simulation, from 2010 to 2100
  time.horizon <-  3
  
  ## Number of runs (i.e. replicas)
  nrun <- 1
  
  ## Initialize model global parameters (equal for all scn)
  spp.distrib.rad <- 20 	# neighborhood radius to determine which species belong to that region (in pixels)
  shrub.colon.rad <- 5 		#
  
  ## Characteristics specific of this scenario
  clim.scn <- "rcp85"
  file.dmnd.harvest <- "DemandHarvest_Bioenergy"
  file.clim.severity <- "ClimaticSeverity_test"
  file.pctg.hot.days <- "PctgHotDays_rcp45"
  file.fire.suppression <- "FireSuppression_CurrExtrem"
    
  ## Spread rate, burn probability parameters, prescribed burns
  rpb <- 1.7
  stochastic.spread <- 1
  pb.upper.th <- 0.8	
  pb.lower.th <- 0.1
  fire.intens.th <- 0.35  # high vs. low intensity fire, SR_noAcc <= fire.intens.th
  pb.target.area <- NA  # if NA, then burnt as 7*pb.convenient.area, otherwise annually burnt pb.fix.area
  pb.convenient.area <- 400 ## should be 15.000 ha
  accum.burnt.area <- rep(pb.convenient.area,7)
  pb.mean <- 1.974
  pb.sd <- 0.683
  pb.fage.th <- 30 ## minimum forest age to apply prescribed burns
  
  
  ## Save all the variables in .r file to be further loaded by landscape.dyn.r
  if(!file.exists(out.path))
    dir.create(file.path(getwd(), out.path), showWarnings = T) 
  dump(ls(), paste0(out.path, "/scn.def.r"))
  
}