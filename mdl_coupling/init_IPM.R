
##should I remove this and only put it in mother R?
#rm(list=ls(all=TRUE))
#gc()
#assign("last.warning", NULL, envir = baseenv())
##

#wd <- "C:/Users/uriso/Desktop/copia_IPM/"
# wd <- getwd()

# IPM functions.
#setwd(wd)
library(tidyverse)
source("./IPM/auxiliary_functions_v9 - Roberto.R")
source("./IPM/IPM_functions_v20_old.R")
library(Rcpp)
sourceCpp("./IPM/IPM_functions_v21_1_Year.cpp")
library(ff)
options(fftempdir = "./")

###########################################
# Read data

# R is control, no change. Other options are "A2" and "B2" (check climate files)
scenario <- "R"

# do we want just tree distribution and abundance or also DBH and demographic dynamics
# note that this can potentially consume a lot of memory

ALL.RESULTS <- F

# final year of simulation. 
# NOTE: if > 2090, it is not compatible with ALL.RESULTS
# NOTE: > 2090 only valid for the CONSTANT CLIMATE SCENARIO (R)

final.year <- 2001

##### KEEP TRACK OF THE SIMULATION NUMBER
##### UPDATE WITH EACH NEW SIMULATION

SIM_NUMBER <- "v39-SP"
sp.tracked <- 1:16

# what do we want to store?
store.original.dbh.dist <- F #remember to change
store.original.abundances <- F #remember to change
store.decadal.dbh.dist <- F
store.decadal.abundances <- F #remember to change

#### colonization parameters
####
max.dist <- 1500

colonization.threshold <- 0.05

# number of saplings to colonize a given plot
# new.saplings <- 100

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

load.workspace <- F

if(!load.workspace){
  
  BA_threshold <- read.csv2("./IPM/BASAL_AREA_THRESHOLD_v3.csv",header = T,stringsAsFactors = F)
  #BA_threshold$perc_95[which(is.na(BA_threshold$perc_95))] <- 15
  BA_threshold$perc_95 <- 15
  
  ##### do we want a log file?
  
  STORE_LOG <- F
  if(STORE_LOG) sink(paste("./IPM/resultados/log_",SIM_NUMBER,".txt",sep=""))
  
  ##### do we want colonization?
  
  COLONIZATION <- T
  
  ##### do we want basal area threshold?
  
  BASAL.AREA.THRESHOLD <- T
  
  map <- read.table("./IPM/MAP_v6.csv",header=T,sep=";",dec=".")

  bcn.map <- read.table(file="./IPM/MAP_BCN_v1.csv",header=T,sep=";",dec=".")
  # cat.map <- read.table(file="./MAP_CAT_v2.csv",header=T,sep=";",dec=".")
  
  map <- subset(map, valid)
  
  # here I trim the area to the province of BCN
  map <- subset(map, ID %in% bcn.map$ID)
  
  # shapefile for plotting
  # select the appropriate one depending on the area
  # either spain, catalunya, or barcelona
  
  # spain <- readOGR(dsn=".",layer="inland.spain_UTM")
  # spain <- readOGR(dsn=".",layer="CAT_UTM_v2")
  spain <- rgdal::readOGR(dsn="./IPM",layer="bcn_UTM_v1") #change?
  
  # necessary for the plots
  spain.df <- fortify(spain)
  ###########################################
  
  ### if the files are already sorted by species, i.e. with a correct "num.sp" field,
  ### we do not need to load the files in any special manner
  
  trees.sorted <- TRUE
  
  if(trees.sorted){
    sp.list <- read.csv2("./IPM/lista_especies_v2.csv",header = T,stringsAsFactors = F)
    newspecies <- unique(sp.list[,3])
    newspecies <- newspecies[trim(newspecies)!=""]
    tesauro <- data.frame(num.sp = c(1:length(newspecies)),name.sp = as.character(unique(sp.list[,3])))
    trees.orig <- read.csv2("./IPM/PIES_MAYORES_IFN3_v10.csv",header = T,stringsAsFactors = F)
    trees.orig$dbh <- as.numeric(trees.orig$dbh)
    trees.orig$factor <- as.numeric(trees.orig$factor)
    
    # plots recorded in the map
    trees.orig <- trees.orig[which(trees.orig$ID %in% map$ID),]
    
    #   # invalid plots: those with presence of plantation, riverside species.
    #   # we do not model these species, and we exclude plots where they are present.
    #   
    #   invalid.id <- unique(append(map$ID[which(map$LAND_USE %in% c(1,2,5))],unique(trees.orig$ID[which(trees.orig$num.sp == -1)])))
    #   
    #   trees.orig <- trees.orig[which(!(trees.orig$ID %in% invalid.id)),]
    
    saplings.orig <- read.csv2("./IPM/SAPLINGS_IFN3_v10.csv",header = T,stringsAsFactors = F)
    saplings.orig <- saplings.orig[which(saplings.orig$num.sp %in% tesauro[,1]),]
    
    saplings.orig <- saplings.orig[which(saplings.orig$ID %in% map$ID),]
    #   saplings.orig <- saplings.orig[which(!(saplings.orig$ID %in% invalid.id)),]
    
    # we only need the plots with presence of saplings. In the original files,
    # plots with seedlings but no saplings are stored.
    
    saplings.orig <- saplings.orig[which(saplings.orig$saplings > 0),]
    
    NUM_SP <- length(newspecies)
    
    ########### ONLY TREES OF SP 5
    #trees.orig <- trees.orig[trees.orig$num.sp == MY.SP,]
    #saplings.orig <- saplings.orig[saplings.orig$num.sp == MY.SP,]
    ###########
    
  }else{
    sp.list <- Readcsv("./IPM/lista_especies_v2")
    newspecies <- levels(unique(sp.list[,3]))
    newspecies <- newspecies[trim(newspecies)!=""]
    tesauro <- data.frame(num.sp = c(1:length(newspecies)),name.sp = as.character(levels(sp.list[,3])))
    NUM_SP <- length(newspecies)
    trees.orig <- load.trees.IFN(3,sp.list)
    saplings.orig <- load.saplings.IFN(3,sp.list)
  }
  
  NUM_PLOTS <- nrow(map)
  
  ###################### some constants
  
  year <- 2000
  
  # IPM Parameters.
  n.time.intervals <- 20
  t.diff <- 10
  initialize.num <- vector("list",NUM_SP)
  
  min.DBH <- 7.5  # in cm.
  n.intervals.mesh <- 1500 ##Discretization DBH
  if (n.intervals.mesh<7) stop("Mesh size must be larger than 8")
  
  #x(cm) is the interval of dbh for every sp.
  x <- x.per.species(min.dbh=min.DBH,n.intervals=n.intervals.mesh)
  y <- x
  x2 <- (x/200)^2
  max.diam <- sapply(1:NUM_SP, function(i) max(x[,i]))
  h <- x[2,]-x[1,]
  nx <- n.intervals.mesh+1
  
  # This matrix is used inside the growth function to save CPU time
  # in calculating the quadrature.
  #### A FASTER WAY?!?!?!?!?!?
  y.minus.x <- array(0,dim=c(dim(x)[1],dim(y)[1],NUM_SP))
  for (k in 1:NUM_SP) {
    for (i in 1:dim(x)[1]) for (j in 1:dim(y)[1]) y.minus.x[i,j,k] <- y[j,k] - x[i,k]
  }
  y.minus.x[y.minus.x<0] <- 0
  
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
  
  ###################### climate
  
  clima <- read.table(file="./IPM/clima/CLIMA_2000.csv",header=T,sep=";",dec=".")
  
  clima.ref <- read.table(file="./IPM/clima/clima_referencia_v5.csv",header=T,sep=";",dec=".")
  
  # double check...
  
  clima <- clima[which(clima$ID %in% map$ID),]
  clima.ref <- clima.ref[which(clima.ref$ID %in% map$ID),]
  
  my.id <- match(map$ID,clima$ID)
  my.id.ref <- match(map$ID,clima.ref$ID)
  
  ###################### 
  
  # trees: ff 3d array
  # saplings: matrix
  # basal area: matrix
  
  adult.trees <- list()
  for(i.sp in 1:length(newspecies)){
    # adult.trees[[newspecies[i.sp]]] <- ff(0,dim = c(NUM_PLOTS,nx))
    adult.trees[[newspecies[i.sp]]] <- matrix(0,NUM_PLOTS,nx)
  }
  ba <- matrix(0,nrow = NUM_PLOTS,ncol = NUM_SP)
  saplings <- matrix(0,nrow = NUM_PLOTS,ncol = NUM_SP)
  
  ################## ALL RESULTS
  ##################
  
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
  ###################### IPM functions coefficients
  
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
  
  load("./IPM/regresiones/colonization_v3")
  conifers <- c(1,4,5,6,7,8,9,10)
  deciduous <- c(2,3)
  quercus <- c(11,12,13,14,15,16)
  
  ###################### other data structures
  
  #load("distancias_bcn_v1")
  load("./IPM/distancias_hasta_2236_v2")
  
  ###############################################################
  ###############################################################
  
  ### first, fill up the trees and saplings distributions
  # for that, we infer a probability density function 
  # from the discrete original data
  orig.adult.trees.file <- paste("./IPM/initial_variables/trees_", "BCN","_", scenario, "_","x_", n.intervals.mesh, ".rdata",sep="")
  orig.ba.file <- paste("./IPM/initial_variables/ba_", "BCN","_", scenario, ".rdata", sep="")
  orig.saplings.file <- paste("./IPM/initial_variables/saplings_", "BCN","_", scenario, ".rdata", sep="")
  
  if (file.exists(orig.adult.trees.file) & file.exists(orig.ba.file) & file.exists(orig.saplings.file)){
    load(orig.adult.trees.file); load(orig.ba.file); load(orig.saplings.file) ##Maybe don't need to load since they are loaded again in IPM_step
    cat("initial variables loaded\n")
  }
  else{
    print(paste(date()," - ",SIM_NUMBER," - ",scenario," - converting IFN data to continuous pdf...",sep=""))
    
    aux.data <- data.frame(ID = sort(rep(map$ID,NUM_SP)),
                           num.sp = rep(1:NUM_SP,NUM_PLOTS),
                           num.plot = sort(rep(1:NUM_PLOTS,NUM_SP)), # just an auxiliary index
                           predicted.basal.area = 0)
    
    trees <- trees.orig %>% group_by(ID,num.sp) %>% summarize(adult.trees = sum(factor))
    aux.data <- merge(aux.data,trees,by.x = c("ID","num.sp"),by.y = c("ID","num.sp"), all.x = T)
    aux.data$adult.trees[is.na(aux.data$adult.trees)] <- 0
    
    aux.data <- merge(aux.data,saplings.orig,by.x = c("ID","num.sp"),by.y = c("ID","num.sp"), all.x = T)
    aux.data$saplings[is.na(aux.data$saplings)] <- 0
    aux.data <- aux.data[,c("ID","num.sp","num.plot","adult.trees","saplings","predicted.basal.area")]
    
    dbh.my.sp <- list()
    for(i.sp in 1:NUM_SP){
      dbh.my.sp[[i.sp]] <- 7.5 + c(0:n.intervals.mesh)*(MAXDBH[i.sp]-7.5)/n.intervals.mesh
    }

    for(i in 1:nrow(aux.data)){
      if(aux.data$adult.trees[i] > 0){
        sp <- aux.data$num.sp[i]
        aux.index <- which(trees.orig$num.sp == aux.data$num.sp[i] & trees.orig$ID == aux.data$ID[i])
        orig.data <- rep(trees.orig$dbh[aux.index],round(trees.orig$factor[aux.index]))
        adult.trees[[sp]][aux.data$num.plot[i],] <- EstimateKernel(orig.data = orig.data,
                                                                   sp = sp,
                                                                   nx = nx,
                                                                   h = h,
                                                                   range.x = c(min(dbh.my.sp[[sp]]),max(dbh.my.sp[[sp]])),
                                                                   bandwidth.range = c(0.03,0.5))
        aux.data$predicted.basal.area[i] <- quadTrapezCpp_1(adult.trees[[sp]][aux.data$num.plot[i],]*x2[,sp],h[sp],nx)*pi
        
      }else{adult.trees[[aux.data$num.sp[i]]][aux.data$num.plot[i],] <- 0}
      if (round(i/10000)*10000==i) print(paste(date(),"- kernel smoothing: running loop",i,"of",nrow(aux.data),sep=" "))
      
    }# for i

    plot.predicted.ba <- aux.data %>% group_by(ID) %>% summarize(ba = sum(predicted.basal.area))
    
    for(i in 1:nrow(aux.data)){
      if(aux.data$predicted.basal.area[i]>0 & aux.data$saplings[i]>0){
        
        my.plot <- aux.data$num.plot[i]
        my.sp <- aux.data$num.sp[i]
        
        q <- c(clima$precip[my.plot],clima$temp[my.plot],plot.predicted.ba$ba[my.plot],clima$anom_pl[my.plot],clima$anom_temp[my.plot])
        
        param.sapl2 <- CoefTimesVarCpp(param=saplings.coef,
                                       z=q,
                                       # s=aux.data$saplings[aux.data$num.plot == my.plot],
                                       spba=aux.data$predicted.basal.area[aux.data$num.plot == my.plot],
                                       type="recruitment")
        
        saplings[my.plot,my.sp] <- ifelse(IPMSaplingsCpp(aux.data$saplings[i],c(param.sapl1[my.sp],param.sapl2[my.sp]))>100,
                                          100,
                                          IPMSaplingsCpp(aux.data$saplings[i],c(param.sapl1[my.sp],param.sapl2[my.sp])))
        if(saplings[my.plot,my.sp] < 0){
          saplings[my.plot,my.sp] <- 0
        }
      }# if basal area and saplings
        if (round(i/10000)*10000==i) print(paste(date(),"- saplings: running loop",i,"of",nrow(aux.data),sep=" "))
    }# for plots
    
    # put back ba values
    for(i.plot in 1:NUM_PLOTS){
      ba[i.plot,] <- aux.data$predicted.basal.area[aux.data$num.plot == i.plot]
    }
    
    # delete auxiliary structures
    
    rm(plot.predicted.ba)
    rm(aux.data)
    rm(trees)
    print(paste(date()," - ",SIM_NUMBER," - ",scenario," - converting IFN data to continuous pdf... done!",sep=""))
    
    save(adult.trees, file=orig.adult.trees.file)
    save(ba, file=orig.ba.file)
    save(saplings, file=orig.saplings.file)
  } #if initial variables files don't exist
  #########################################################
  #########################################################

  if(store.original.dbh.dist){
    save(list=ls(all=T),file=paste(scenario,"_",SIM_NUMBER,".RData",sep=""))
    for(i.sp in 1:length(newspecies)){
      temp <- adult.trees[[i.sp]]
      # ffsave(temp,file = paste("./resultados/trees_array_original_sp_",i.sp,"_",SIM_NUMBER,"_",scenario,sep=""))
      save(temp,file = paste("./IPM/resultados/trees_array_original_sp_",i.sp,"_",SIM_NUMBER,"_",scenario,sep=""))
    }
    rm(temp)  
  }
  ### store original results
  
  if(store.original.abundances){
    SaveAbundance(map = map,
                  trees = adult.trees,
                  saplings = saplings,
                  tesauro = tesauro,
                  year = "IFN3",
                  scenario = scenario,
                  SIM_NUMBER = SIM_NUMBER,
                  plot.distribution = F,
                  plot.abundance = F,
                  plot.richness = F,
                  plot.ba = F,
                  spain.df,
                  path="./IPM/resultados/")
  }
  
}else{
  load(file=paste("./IPM/",scenario,"_",SIM_NUMBER,".RData",sep=""))
  adult.trees <- list()
  for(i.sp in 1:16){
    ffload(file = paste("./IPM/resultados/trees_array_original_sp_",i.sp,"_",SIM_NUMBER,"_",scenario,sep=""),overwrite = T)
    adult.trees[[i.sp]] <- temp
    rm(temp)
  }
}
adult.trees.file <- paste("./mdl_coupling/dyn_var_IPM/trees_", "BCN","_", scenario, "_","x_", n.intervals.mesh, ".rdata",sep="")
ba.file <- paste("./mdl_coupling/dyn_var_IPM/ba_", "BCN","_", scenario, ".rdata", sep="")
saplings.file <- paste("./mdl_coupling/dyn_var_IPM/saplings_", "BCN","_", scenario, ".rdata", sep="")