library(sp)
library(raster)  
library(RANN)  # for nn2()
library(tidyverse)
source("mdl/update.clim.r")
source("mdl/update.interface.r")
source("mdl/prob.igni.r")
source("mdl/growth.r")
source("mdl/drought.r")
source("mdl/cohort.establish.r")
source("mdl/afforestation.r")
source("mdl/forest.mgmt.r")
source("mdl/fire.regime.r")
source("mdl/post.fire.r")
source("mdl/auxiliars.r")

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
#load("inputlyrs/rdata/harvest.rdata")
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


## Climatic severity and pctg hot days tables
clim.severity <- read.table(paste0("inputfiles/", file.clim.severity, ".txt"), header=T)


## Build the baseline time sequence and the time sequence of the processes (shared for all runs). 
## 1. Climate change, 2. Land-cover changes, 3. Forest management
## 4. Wildfires, 5. Prescribed burns, 6. Drought, 7. Post-fire regeneration,
## 8. Cohort establihsment, 9. Afforestation, 10. Growth
time.seq <- seq(1, time.horizon, 1)

if(time.horizon==1){
  clim.schedule <- 1
  } else{
  clim.schedule <- seq(1, time.horizon-1, clim.step)}
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
track.fire <-  data.frame(run=NA, year=NA, swc=NA, clim.sever=NA, fire.id=NA, fst=NA, 
                          wind=NA, atarget=NA, aburnt.highintens=NA, 
                          aburnt.lowintens=NA, asupp.fuel=NA, asupp.sprd=NA)
track.pb <-  data.frame(run=NA, year=NA, clim.sever=NA, fire.id=NA, 
                        wind=NA, atarget=NA, aburnt.lowintens=NA)
track.drougth <- data.frame(run=NA, year=NA, spp=NA, ha=NA)
track.cohort <- data.frame(run=NA, year=NA, spp.out=NA, Var2=NA, Freq=NA)
track.post.fire <- data.frame(run=NA, year=NA, spp.out=NA, Var2=NA, Freq=NA)
track.afforest <- data.frame(run=NA, year=NA, Var1=NA, Freq=NA)
track.land <- data.frame(run=NA, year=NA, spp=NA, area=NA, vol=NA, volbark=NA, carbon=NA)