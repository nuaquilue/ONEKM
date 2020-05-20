source(paste0(work.path,"/mdl_interface/select.valid.plots.r"))
cat("Initializing IPM parameters", "\n")

##study area
study.area <- "BCN"
##load valid plots map
if(file.exists(paste0("/mdl_interface/IPM_",study.area, "_map.rdata"))){
	load(paste0("/mdl_interface/IPM_",study.area, "_map.rdata"))
	cat("IPM valid plots map loaded\n")
} else{
	cat("creating IPM valid plots map\n")
	select.valid.plots(study.area)
	load((paste0(work.path,"/mdl_interface/IPM_",study.area, "_map.rdata"))
}
NUM_PLOTS <- nrow(map)
# spain <- readOGR(dsn=".",layer="inland.spain_UTM")
# spain <- readOGR(dsn=".",layer="CAT_UTM_v2")
spain <- rgdal::readOGR(dsn="./IPM",layer="bcn_UTM_v1")
## necessary for the plots
spain.df <- fortify(spain)

## what do we want to store?
save.IPM.variables <- F
store.original.dbh.dist <- F 
store.original.abundances <- F 
store.decadal.dbh.dist <- F
store.decadal.abundances <- F
## do we want just tree distribution and abundance or also DBH and demographic dynamics
## note that this can potentially consume a lot of memory
ALL.RESULTS <- F

##initial year, variables are initialised at 2000 so unless initialisation is modified don't change
year <- 2000

##species to study
sp.list <- read.csv2("./IPM/lista_especies_v2.csv",header = T,stringsAsFactors = F)
newspecies <- unique(sp.list[,3])
newspecies <- newspecies[trim(newspecies)!=""]
tesauro <- data.frame(num.sp = c(1:length(newspecies)),name.sp = as.character(unique(sp.list[,3])))
NUM_SP <- length(newspecies)
sp.tracked <- 1:16
conifers <- c(1,4,5,6,7,8,9,10)
deciduous <- c(2,3)
quercus <- c(11,12,13,14,15,16)

## adult trees discretization parameters
min.DBH <- 7.5  # in cm.
n.intervals.mesh <- 4000 ##Discretization DBH
x <- x.per.species(min.dbh=min.DBH,n.intervals=n.intervals.mesh)
y <- x
x2 <- (x/200)^2
max.diam <- sapply(1:NUM_SP, function(i) max(x[,i]))
h <- x[2,]-x[1,]
nx <- n.intervals.mesh+1
MAXDBH <- numeric(16)
MAXDBH[1] <- 155.54    # conifers.
MAXDBH[2] <- 290.73    # Deciduous.
MAXDBH[3] <- 330   # Fagus sylvatica.
MAXDBH[4] <- 220   # Juniperus thurifera.
MAXDBH[5] <- 160.6   # Pinus halepensis.
MAXDBH[6] <- 162.8   # Pinus nigra.
MAXDBH[7] <- 198   # Pinus pinaster.
MAXDBH[8] <- 150.7   # Pinus pinea.
MAXDBH[9] <- 163.9   # Pinus sylvestris.
MAXDBH[10] <- 141.9  # Pinus uncinata.
MAXDBH[11] <- 198  # Quercus faginea.
MAXDBH[12] <- 167.2  # Quercus ilex.
MAXDBH[13] <- 189.2  # Quercus pyrenaica.
MAXDBH[14] <- 313.5  # Quercus robur/petraea.
MAXDBH[15] <- 161.7  # Quercus suber
MAXDBH[16] <- 114.84   # Sclerophyllous

##IPM parameters
survival.coef <- load.survival("survival_v8",newspecies)
growth.coef <- load.growth("growth_v6",newspecies)
ingrowth.coef <- load.ingrowth("ingrowth_v10",newspecies)
saplings.coef <- load.saplings("recruitment_regression_v15",newspecies)
param.survival1 <- survival.coef$log.dbh
param.growth1 <- growth.coef$log.dbh
param.growth3 <- growth.coef$intercept.variance
param.growth4 <- growth.coef$slope.variance
param.sapl1 <- saplings.coef$binom
param.ingrowth1 <- ingrowth.coef$lambda

## Colonization parameters
COLONIZATION <- T ##do we want colonization?
max.dist <- 1500 ##meters
colonization.threshold <- 0.05
BASAL.AREA.THRESHOLD <- T
BA_threshold <- read.csv2("./IPM/BASAL_AREA_THRESHOLD_v3.csv",header = T,stringsAsFactors = F)
BA_threshold$perc_95 <- 15
## number of saplings to colonize a given plot
new.saplings <- c(297.9471, 
                  353.4785,
                  264.4421,
                  339.7690,
                  485.3332,
                  324.3320,
                  485.3627,
                  322.2076,
                  322.2275,
                  254.6479,
                  279.6451,
                  516.2391,
                  605.3193,
                  309.4185,
                  242.8038,
                  502.7275)
load("./IPM/regresiones/colonization_v3")
load("./IPM/distancias_hasta_2236_v2") 

##hard coded fire.regeneration, needs to be changed so it's a file
fire.regeneration <- c(T,T,T,F,T,F,T,F,F,F,T,T,T,T,T,T)

orig.adult.trees.file <- paste("./IPM/initial_variables/trees_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
orig.ba.file <- paste("./IPM/initial_variables/ba_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
orig.saplings.file <- paste("./IPM/initial_variables/saplings_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
orig.plots.age.file <- paste("./IPM/initial_variables/orig_plots_age", study.area, ".rdata",sep="")