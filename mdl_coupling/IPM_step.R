
IPM.step <- function(iyear){

  ##########definition of dyn var files#######################
  #should put this in initialisation r
  adult.trees.file <- paste("./dyn_var/trees_", "BCN","_", scenario, "_","x_", n.intervals.mesh, ".rdata",sep="")
  ba.file <- paste("./dyn_var/ba_", "BCN","_", scenario, ".rdata", sep="")
  saplings.file <- paste("./dyn_var/saplings_", "BCN","_", scenario, ".rdata", sep="")
  #############################################################
  
  if (file.exists(adult.trees.file) & file.exists(ba.file) & file.exists(saplings.file)){
    load(adult.trees.file); load(ba.file); load(saplings.file) ##make sure that it updates values
  }
  else{
    load(orig.adult.trees.file); load(orig.ba.file); load(orig.saplings.file)
  }
  
  print(paste(date(),"- starting IPM - year",iyear,sep=" "))
  
  ### Read climate data if necessary
  if(iyear != 2000){
    if(scenario != "R"){
      if(iyear != 2010){
        if(scenario == "A2"){
          clima <- read.table(file=paste("./clima/MEAN_A2_",iyear,".csv",sep=""),header=T,sep=";",dec=".")
        }else{
          clima <- read.table(file=paste("./clima/MEAN_B2_",iyear,".csv",sep=""),header=T,sep=";",dec=".")
        }
      }else{
        clima <- read.table(file="./clima/MEAN_2010.csv",header=T,sep=";",dec=".")
      }
    }
  }
  my.id <- match(map$ID,clima$ID)
  
  if(ALL.RESULTS){
    # indexes for storing results
    ###################################################################################
    ###################################################################################
    my.new.adults <- which(names(results[[1]]) == paste("new_adults_",iyear,sep=""))
    my.adults <- which(names(results[[1]]) == paste("adults_",iyear,sep=""))
    my.deaths <- which(names(results[[1]]) == paste("deaths_",iyear,sep=""))
    my.saplings <- which(names(results[[1]]) == paste("saplings_",iyear,sep=""))
    my.ba <- which(names(results[[1]]) == paste("basal_area_",iyear,sep=""))
    
    previous.adults <- which(names(results[[1]]) == paste("adults_",(iyear - 10),sep=""))
    ###################################################################################
    ###################################################################################
  }

  # IPM loop        
  for (i in 1:NUM_PLOTS) {
    
    if (round(i/100)*100==i) print(paste(date(),"- IPM: running loop",i,"of",NUM_PLOTS," - year",iyear,sep=" "))
    
    sapl <- saplings[i,]
    
    # are there adult trees or saplings?
    if (sum(ba[i,])>0 | sum(sapl)>0){
      
      temper <- clima$temp[my.id[i]]
      rain <- clima$precip[my.id[i]]
      anom.temper <- clima$anom_temp[my.id[i]]  
      anom.rain <- clima$anom_pl[my.id[i]]
      
      q <- c(rain,temper,sum(ba[i,]),anom.rain,anom.temper)
      
      ########################
      ########################
      
      # if adult trees, apply IPM
      if(sum(ba[i,])>0){
        
        sum.ntrees <- numeric(length(newspecies))
        
        for(i.sp in 1:length(newspecies)){
          # this is not the number of trees, it's just for knowing if there are trees of that sp
          sum.ntrees[i.sp] <- sum(adult.trees[[i.sp]][i,])
          #apply(X = adult.trees[[i.sp]][i,],MARGIN = 1,FUN = sum) 
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
                                                     nx=nx,
                                                     y_minus_x=y.minus.x[,,j])          
            adult.trees[[j]][i,] <- dummy
            p.sapl <- c(param.sapl1[j],param.sapl2[j])
            
            # careful with the saplings prediction. If not controlled, it can both explode and drop below 0
            saplings[i,j] <- ifelse(IPMSaplingsCpp(sapl[j],p.sapl)>100,100,IPMSaplingsCpp(sapl[j],p.sapl))
            if(saplings[i,j] < 0){
              saplings[i,j] <- 0
            }
            if(ALL.RESULTS){
              ###################################################################################
              ###################################################################################
              
              # adults and saplings are computed the same way in each step, because
              # of the data conversion to pdfs
              
              results[[j]][i,my.adults] <- quadTrapezCpp_1(dummy,h[j],nx)
              results[[j]][i,my.saplings] <- ifelse(is.null(saplings[i,j]),
                                                    0,
                                                    saplings[i,j]*(10000/(pi*25)))
              
              # deaths, on the other hand, are different for the first step
              if(iyear == 2000){
                results[[j]][i,my.deaths] <- ifelse(sum(trees.orig$factor[trees.orig$num.sp == j & trees.orig$ID == map$ID[i]]) - results[[j]]$adults_2000[i] > 0,
                                                    sum(trees.orig$factor[trees.orig$num.sp == j & trees.orig$ID == map$ID[i]]) - results[[j]]$adults_2000[i],
                                                    0)
              }else{
                results[[j]][i,my.deaths] <- ifelse(results[[j]][i,previous.adults] - results[[j]][i,my.adults] > 0,
                                                    results[[j]][i,previous.adults] - results[[j]][i,my.adults],
                                                    0)
              }
              ###################################################################################
              ###################################################################################
            }
          }#sum(ntrees[[j]])
        }# for j in NUM_SP
      }# if(sum(unlist(trees[[i]]$ba))>0)
      
      ########################
      ########################
      
      # if saplings on the previous step, apply ingrowth IPM
      if(sum(sapl)>0){
        
        param.ingrowth2 <-  CoefTimesVarCpp(param=ingrowth.coef,z=q,s=sapl,type="ingrowth") 
        
        order <- sample(x = 1:NUM_SP,size = NUM_SP,replace = F)
        for(j in order){
          
        if(sapl[j] > 0 & (!BASAL.AREA.THRESHOLD | (sum(ba[i,]) < BA_threshold$perc_95[j]))){

          dummy <- IPMIngrowthIdentityCpp(y[,j],c(param.ingrowth1[j],param.ingrowth2[j]))
          if(ALL.RESULTS){
            ###################################################################################
            results[[j]][i,my.new.adults] <- quadTrapezCpp_1(dummy,h[j],nx)
            ###################################################################################
          }
          adult.trees[[j]][i,] <- adult.trees[[j]][i,] + dummy
          if(ALL.RESULTS){
            ###################################################################################
            results[[j]][i,my.adults] <- quadTrapezCpp_1(adult.trees[[j]][i,],h[j],nx)
            ###################################################################################
          }
          }# if basal area threshold is met for j-th species and there are saplings
        } # for all j-th species
      } # if there are saplings of any sp
      
      ########################
      ########################
      
      # if any of the two conditions is met (saplings or adults), update basal area of the plot
      order <- sample(x = 1:NUM_SP,size = NUM_SP,replace = F)
      for (j in order) {
        dummy <- adult.trees[[j]][i,]
        num.trees <- quadTrapezCpp_1(dummy,h[j],nx)
        if(num.trees < 0.1){
          adult.trees[[j]][i,] <- 0
          ba[i,j] <- 0
        }else{
          ba[i,j] <- quadTrapezCpp_1(dummy*x2[,j],h[j],nx)*pi
        }
        if(ALL.RESULTS){
          ###################################################################################
          results[[j]][i,my.ba] <- ba[i,j] #ifelse(is.null(trees[[i]]$ba[[j]]),0,trees[[i]]$ba[[j]])
          ###################################################################################
        }
      }# for j in NUM_SP
    }# if saplings or adults
    
  }# for i in plots

  #########################################################
  #########################################################
  
  if(COLONIZATION){
    
    # keep a dataframe with the basal area for potential colonization
    temp.suitable <- data.frame(ID = map$ID,
                                saplings = integer(nrow(map)),
                                plot_basal_area = integer(nrow(map)),
                                #neigh = integer(nrow(map)),
                                neigh_ba = integer(nrow(map)),
                                suitable = logical(nrow(map)))
    temp.suitable$plot_basal_area <- rowSums(ba) #sapply(X=trees,FUN=function(x) sum(unlist(x$ba)))
    # default 
    temp.suitable$suitable <- T
    order <- sample(x = 1:NUM_SP,size = NUM_SP,replace = F)
    for(i in order){
      
      # for every species:
      # - select suitable plots
      # - load regression coefficients
      # - predict new number of saplings
      # - append new saplings to the saplings file
      
      # suitable plots are these without presence,
      # closer than 2236 to a plot with presence,
      # with suitable land use,
      # and total basal area lower than the specific threshold
      
      print(paste(date()," - ",SIM_NUMBER," - ",scenario," - year ",iyear," - colonization for sp ",i,"...",sep=""))
      suitable <- GetSuitablePlots(map=map,
                                   temp.suitable=temp.suitable,
                                   adult.trees=adult.trees,
                                   ba = ba,
                                   #trees.index=trees.index,
                                   distancias=distancias,
                                   sp=tesauro[i,1],
                                   max.dist = max.dist,
                                   verbose=T)   
      
      #apply colonization
      if(!is.null(nrow(suitable))){
        if(nrow(suitable)>0){
          if(i %in% conifers) my.model <- colonization.glm[[1]]
          else if(i %in% quercus) my.model <- colonization.glm[[2]]
          else if(i %in% deciduous) my.model <- colonization.glm[[3]]
          
          suitable$colonized <- predict(my.model,newdata = suitable,type = "response")
          suitable$colonized <- ifelse(suitable$colonized > colonization.threshold,1,0)
          
          k <- match(suitable$ID,map$ID)
          #add new saplings to saplings list
          for(j in 1:nrow(suitable)){
            # saplings[k[j],tesauro[i,1]] <- suitable$saplings[j] + saplings[k[j],tesauro[i,1]]
            saplings[k[j],tesauro[i,1]] <- new.saplings[tesauro[i,1]] + saplings[k[j],tesauro[i,1]]
          }
        }
      }# if there are suitable plots
    }#for sp.
  }
  #########################################################
  #########################################################
  print(paste(date()," - ",SIM_NUMBER," - ",scenario," - year ",iyear," - saving results...",sep=""))
  
  ### store results
  if(store.decadal.abundances){
    SaveAbundance(map = map,
                  trees = adult.trees,
                  saplings = saplings,
                  tesauro = tesauro,
                  year = iyear,
                  scenario = scenario,
                  SIM_NUMBER = SIM_NUMBER,
                  plot.distribution = T,
                  plot.abundance = T,
                  plot.richness = F,
                  plot.ba = F,
                  spain.df,
                  path="./resultados/")
  }
  if(store.decadal.dbh.dist){
    save(list=ls(all=T),file=paste(SIM_NUMBER,"_",scenario,"_",iyear,".RData",sep=""))
    for(i.sp in 1:length(newspecies)){
      temp <- adult.trees[[i.sp]]
      ffsave(temp,file = paste("./resultados/trees_array_year_",iyear,"_sp_",i.sp,"_",SIM_NUMBER,"_",scenario,sep=""))
    }
    rm(temp)  
  }
  save(adult.trees, file=adult.trees.file)
  save(ba, file=ba.file)
  save(saplings, file=saplings.file)

}#for iyear