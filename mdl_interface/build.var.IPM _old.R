build.var.IPM <- function(n.intervals.mesh, study.area){
    
    if (n.intervals.mesh<7) stop("Mesh size must be larger than 8")

    ##Reads species to be studied
    sp.list <- read.csv2("./IPM/lista_especies_v2.csv",header = T,stringsAsFactors = F)
    newspecies <- unique(sp.list[,3])
    newspecies <- newspecies[trim(newspecies)!=""]
    NUM_SP <- length(newspecies)

    ##Discretisation parameters
    min.DBH <- 7.5  # in cm.
    x <- x.per.species(min.dbh=min.DBH,n.intervals=n.intervals.mesh)
    x2 <- (x/200)^2
    h <- x[2,]-x[1,]
    nx <- n.intervals.mesh+1

    ##Reads IFN data of plots to study
    load(paste0(work.path,"/mdl_interface/map.rdata"))
    NUM_PLOTS <- nrow(map)
    tesauro <- data.frame(num.sp = c(1:length(newspecies)),name.sp = as.character(unique(sp.list[,3])))
    trees.orig <- read.csv2("./IPM/PIES_MAYORES_IFN3_v10.csv",header = T,stringsAsFactors = F)
    trees.orig$dbh <- as.numeric(trees.orig$dbh)
    trees.orig$factor <- as.numeric(trees.orig$factor)
    trees.orig <- trees.orig[which(trees.orig$ID %in% map$ID),]

    saplings.orig <- read.csv2("./IPM/SAPLINGS_IFN3_v10.csv",header = T,stringsAsFactors = F)
    saplings.orig <- saplings.orig[which(saplings.orig$num.sp %in% tesauro[,1]),]
    saplings.orig <- saplings.orig[which(saplings.orig$ID %in% map$ID),]
    saplings.orig <- saplings.orig[which(saplings.orig$saplings > 0),]        

    ##INITIALISATION OF IPM VARIABLES
    ## trees (adult.trees): list of nºspecies matrices -> rows:plots and colums:dbh
    ## saplings: matrix -> rows:plots and columns:species
    ## basal area (ba): matrix -> rows:plots and columns:species
    ## Matrix plots indexes are the same as in map where they can be associated to plot ID

    orig.adult.trees.file <- paste("./IPM/initial_variables/trees_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
    orig.ba.file <- paste("./IPM/initial_variables/ba_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
    orig.saplings.file <- paste("./IPM/initial_variables/saplings_", study.area,"_", n.intervals.mesh, ".rdata",sep="")

    if (file.exists(orig.adult.trees.file) & file.exists(orig.ba.file) & file.exists(orig.saplings.file)){
      load(orig.adult.trees.file); load(orig.ba.file); load(orig.saplings.file) ##Maybe don't need to load since they are loaded again in IPM_step
      cat("initial variables loaded\n")
    
    } else{

      adult.trees <- list()
      for(i.sp in 1:length(newspecies)){
      # adult.trees[[newspecies[i.sp]]] <- ff(0,dim = c(NUM_PLOTS,nx))
      adult.trees[[newspecies[i.sp]]] <- matrix(0,NUM_PLOTS,nx)
      }
      ba <- matrix(0,nrow = NUM_PLOTS,ncol = NUM_SP)
      saplings <- matrix(0,nrow = NUM_PLOTS,ncol = NUM_SP)

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

    }
    
    ##Calculate empty plots that LFT assigns as forested plots
    empty.plots.forest <- c()
    for (i in 1:nrow(ba)){
      if(sum(ba[i,]==0 & map$spp[i])<14){
        empty.plots.forest <- c(empty.plots.forest, i)
      }
    }

    ##Builds forest in the forested empty cells by averaging the trees of its 8 closest neighbours
    aux.adult.trees <- list()
    for(i.sp in 1:length(newspecies)){
      aux.adult.trees[[newspecies[i.sp]]] <- matrix(0,length(empty.plots.forest),nx)
    }
    aux.saplings <- matrix(0,nrow = length(empty.plots.forest),ncol = NUM_SP)
    aux.ba <- matrix(0,nrow = length(empty.plots.forest),ncol = NUM_SP)

    default.neigh <- c(289L, 290L, 291L, -289L, -290L, -291L, -1L, 1L)
    for (i in 1:length(empty.plots.forest)){
      plot.neigh <- map$Medfire.id + default.neigh 
      plot.neigh <- which(map$Medfire.id %in% plot.neigh)
      for(j in 1:length(newspecies)){
        aux.adult.trees[[j]][i,] <- rowMeans(adult.trees[[j]][plot.neigh,])
        aux.saplings[i,j] <- mean(saplings[plot.neigh,j])
        aux.ba[i,j] <- quadTrapezCpp_1(aux.adult.trees[[j]][i,]*x2[,j],h[j],nx)*pi
      }
    }

    for(j in 1:length(newspecies)){
      adult.trees[[j]][empty.plots.forest,] <- aux.adult.trees[[j]]
      saplings[empty.plots.forest,j] <- aux.saplings[,j]
      ba[empty.plots.forest,j] <- aux.ba[,j]
    }



    save(adult.trees, file=orig.adult.trees.file)
    save(ba, file=orig.ba.file)
    save(saplings, file=orig.saplings.file)

    print("new initial variables saved")
}