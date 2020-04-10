rm(list=ls())

###Change to local directory###
# setwd("c:/work/MEDMOD/SpatialModelsR/MEDFIRE")  #Nú HP
# setwd("d:/MEDMOD/SpatialModelsR/MEDFIRE")   #CTFC
setwd("C:/Users/uriso/Desktop/ONEKM") #Oriol laptop
## Load required packages and functions 
library(tidyverse)
library(sp)
library(raster)  
library(RANN)  # for nn2()

map.csv <- "MAP_BCN_v1.csv"

# set the scenario
source("mdl/define.scenario.r")
scn.name <- "ONEKM-v2"
define.scenario(scn.name)

#run coupled
source("mdl_coupling/run_coupled.r")
system.time(IPM.Medfire.mdl(scn.name))

# run the model
source("mdl/land.dyn.mdl.r")  
system.time(land.dyn.mdl(scn.name))

# Initialization functions
source("mdl/read.static.vars.r")
source("mdl/read.state.vars.r")
source("mdl/read.climatic.vars.r")
#work.path <- "C:/WORK/MEDMOD/SpatialModels/MEDFIRE_II"
#work.path <- "C:/WORK/MEDMOD/SpatialModelsR/MEDFIRE"
work.path <- "C:/Users/uriso/Desktop/ONEKM"
#adapt.climatic.vars(work.path)
# Create .Rdata with static variables of the model, only run once for all scenarios!
read.static.vars(work.path)
# Create .Rdata with initial values of variables of the model, used at each replicate of any scn.
read.state.vars(work.path)
# Create a data frame per climatic scenario and decade with climatic variables (temp and precip) of a specific model for CAT 
read.climatic.vars2(work.path, "SMHI-RCA4_MOHC-HadGEM2-ES")


## Save interfaces
source("mdl/update.interface.r")
load("inputlyrs/rdata/land.rdata")
interface <- update.interface(land)
save(interface, file="inputlyrs/rdata/interface.rdata")