rm(list=ls())
setwd("c:/work/MEDMOD/SpatialModelsR/MEDFIRE")  #Nú HP
# setwd("d:/MEDMOD/SpatialModelsR/MEDFIRE")   #CTFC

# set the scenario
source("mdl/define.scenario.r")
scn.name <- "TestOnekm"
define.scenario(scn.name)
# run the model
source("mdl/land.dyn.mdl.r")  
system.time(land.dyn.mdl(scn.name))

# Initialization functions
source("mdl/read.static.vars.r")
source("mdl/read.state.vars.r")
source("mdl/read.climatic.vars.r")
work.path <- "C:/WORK/MEDMOD/SpatialModels/MEDFIRE_II"
work.path <- "C:/WORK/MEDMOD/SpatialModelsR/MEDFIRE"

# Create .Rdata with static variables of the model, only run once for all scenarios!
read.static.vars(work.path)
# Create .Rdata with initial values of variables of the model, used at each replicate of any scn.
read.state.vars(work.path)
# Create a data frame per climatic scenario and decade with climatic variables (temp, precip and rad) for CAT
read.climatic.vars(work.path)

## Save interfaces
source("mdl/update.interface.r")
load("inputlyrs/rdata/land.rdata")
interface <- update.interface(land)
save(interface, file="inputlyrs/rdata/interface.rdata")