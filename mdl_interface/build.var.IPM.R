build.var.IPM <- function(n.intervals.mesh, study.area){

    load(map.adress)#to do
    sp.list <- read.csv2("./IPM/lista_especies_v2.csv",header = T,stringsAsFactors = F)
    newspecies <- unique(sp.list[,3])
    newspecies <- newspecies[trim(newspecies)!=""]
    tesauro <- data.frame(num.sp = c(1:length(newspecies)),name.sp = as.character(unique(sp.list[,3])))
    trees.orig <- read.csv2("./IPM/PIES_MAYORES_IFN3_v10.csv",header = T,stringsAsFactors = F)
    trees.orig$dbh <- as.numeric(trees.orig$dbh)
    trees.orig$factor <- as.numeric(trees.orig$factor)
    trees.orig <- trees.orig[which(trees.orig$ID %in% map$ID),]

    saplings.orig <- read.csv2("./IPM/SAPLINGS_IFN3_v10.csv",header = T,stringsAsFactors = F)
    saplings.orig <- saplings.orig[which(saplings.orig$num.sp %in% tesauro[,1]),]
    saplings.orig <- saplings.orig[which(saplings.orig$ID %in% map$ID),]
    saplings.orig <- saplings.orig[which(saplings.orig$saplings > 0),]

    NUM_SP <- length(newspecies)
        
    if (n.intervals.mesh<7) stop("Mesh size must be larger than 8")

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

    orig.adult.trees.file <- paste("./IPM/initial_variables/trees_", "BCN","_", scenario, "_","x_", n.intervals.mesh, ".rdata",sep="")
    orig.ba.file <- paste("./IPM/initial_variables/ba_", "BCN","_", scenario, ".rdata", sep="")
    orig.saplings.file <- paste("./IPM/initial_variables/saplings_", "BCN","_", scenario, ".rdata", sep="")

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

    construct.forests()#to do

    save(adult.trees, file=orig.adult.trees.file)
    save(ba, file=orig.ba.file)
    save(saplings, file=orig.saplings.file)

}