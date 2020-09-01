###
### Auxiliary functions for dispersal scripts
###
### v1: JDGC - 03/2012
###
### basal_area(pmfile,idparcela,IFN,sp)
### surrounding_basal_area(pmfile,mapfile,idparcela,mindist,maxdist,IFN,sp)
### surrounding_stands(pmfile,mapfile,idparcela,mindist,maxdist,IFN,sp)
### colonization_cat_saplings25(stands,idparcela)
### colonization_saplings75(stands,idparcela)
### colonization_pm(stands,idparcela)
### UTMDistance(stand1,stand2)
### UTMDistance2(stand1,X,Y)
### IsFront(stands,idparcela,IFN)

#### v9 - reprogram SaveAbundance and GetSuitablePlots
#### to accept adult.trees[[NUM_SP]][NUM_PLOTS,nx]
#### where each [NUM_PLOTS,nx] matrix is an ff object

#library(plyr)
library(ggplot2)
library(MASS)
library(emdbook) ##remember to allow for working version
library(foreach)
# library(gregmisc)
library(raster)
library(nlme)
library(bbmle)
library(reshape2)
library(grid)
library(gridExtra)
library(rgdal)
library(mclust)
library(KernSmooth)
library(dplyr)

# library(doMC)
# registerDoMC(6)

#########

Writecsv <- function(myfile,name){
  write.table(myfile,file=paste(name,".csv",sep=""),sep=";",dec=".",append=FALSE,row.names=FALSE,col.names=TRUE)
}

#######################################
#######################################

WritecsvFormat <- function(myfile,name){
  write.table(format(myfile,scientific=FALSE,digits=6),file=paste(name,".csv",sep=""),quote=FALSE,sep=";",dec=".",append=FALSE,row.names=FALSE,col.names=TRUE)
}

#######################################
#######################################

Readcsv <- function(myfile){
  read.csv2(paste(myfile,".csv",sep=""),header=TRUE,dec=".",sep=";",comment.char="")
}

#######################################
#######################################

GetCoord <- function(x,y){
  return (paste(as.character(x),as.character(y),sep=""))
}


#######################################
#######################################

GetIndex <- function(list){
  index <- matrix(nrow=length(list),ncol=2)
  index[,1] <- c(1:length(list))
  test <- unlist(list[])
  test <- test[which(names(test)=="ID")]
  index[,2] <- test
  rm(test)
  gc()
  return(index)
}

#######################################
#######################################

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#######################################
#######################################

EstimateKernel <- function(orig.data,
                           sp,
                           nx,
                           h,
                           range.x,
                           bandwidth.range = c(0.01,0.5)){
  bandwidth <- try(dpik(orig.data),silent=T)
  # if (class(bandwidth)=="try-error") {
  #   bandwidth <- try(bw.SJ(x=orig.data,nb=nx,method="dpi"),silent=T)
  #   if (class(bandwidth)=="try-error") {
  #     best.bw <- 0
  #     best.value <- 10000
  #     for(i.bw in seq(bandwidth.range[1],bandwidth.range[2],0.01)){
  #       my.kernel <- suppressWarnings(bkde(x=orig.data,bandwidth=i.bw,gridsize=nx,range.x=range.x))
  #       my.kernel <- my.kernel$y*length(orig.data)
  #       estimated.count <- quad.trapez(my.kernel,h[sp],nx,1)
  #       if(estimated.count - length(orig.data)<best.value){
  #         best.value <- estimated.count - length(orig.data)
  #         best.bw <- i.bw
  #       }# if better
  #     }# for bandwidth values
  #     bandwidth <- best.bw
  #   }# if bw.SJ fails
  # }#if dpik fails 
  
# Modification 06/02/2020.
  bandwidth <- .5
  
  # multiply kernel by number of trees, as to get the integral right
  out <- bkde(x=orig.data,bandwidth=bandwidth,gridsize=nx,range.x=range.x)$y*length(orig.data)
  return(out)
}

#######################################
#######################################

EstimateKernel_iterative <- function(orig.data,sp,nx,range.x,bandwidth.range = c(0.01,0.5)){

      best.bw <- 0
      best.value <- 10000
      for(i.bw in seq(bandwidth.range[1],bandwidth.range[2],0.01)){
        my.kernel <- suppressWarnings(bkde(x=orig.data,bandwidth=i.bw,gridsize=nx,range.x=range.x))
        my.kernel <- my.kernel$y*length(orig.data)
        estimated.count <- quad.trapez(my.kernel,h[sp],nx,1)
        if((estimated.count - length(orig.data))<best.value){
          best.value <- estimated.count - length(orig.data)
          best.bw <- i.bw
        }# if better
      }# for bandwidth values
      bandwidth <- best.bw
    
  # multiply kernel by number of trees, as to get the integral right
  out <- bkde(x=orig.data,bandwidth=bandwidth,gridsize=nx,range.x=range.x)$y*length(orig.data)
  return(out)
}

#######################################
#######################################

quad.trapez <- function(q,h,nx,dim.int=2) {
  if (dim.int==1) int.q <- sum(q) - .5*(q[1]+q[nx])
  else int.q <- colSums(q) - .5*(q[1,]+q[nx,])
  return(h*int.q)
}

interval.quad.trapez <- function(q,h,dim.int=1) {
  if (dim.int==1) int.q <- sum(q) - .5*(q[1]+q[2])
  return(h*int.q)
}

helper.get.adults <- function(trees.list,num.sp){
  ifelse(length(trees.list$num[[num.sp]])<=1,0,quad.trapez(trees.list$num[[num.sp]],h[num.sp],nx,1))
}

helper.get.spba <- function(trees.list,num.sp){
  ifelse(is.null(trees.list$ba[[num.sp]]),0,trees.list$ba[[num.sp]])
}

helper.get.saplings <- function(sap.list,num.sp){
  ifelse(is.null(sap.list$num[[num.sp]]),0,sap.list$num[[num.sp]])
}

helper.get.ba <- function(trees.list){
  sum(unlist(trees.list$ba))
}

#######################################
#######################################

SaveDBH <- function(map,
                    trees,
                    type="list",
                    nx,
                    sp,
                    max.size,
                    sim.number = "0",
                    scenario = "R",
                    year = 0,
                    path = "/home/david/CREAF/articulo/IPM/resultados/",
                    plot.index = 0,
                    write.to.file=F){
  
  # this is the array we will save
  dbh.distribution <- numeric(nx)
  
  # which are the intervals for our particular species
  dbh.my.sp <- 7.5 + c(0:nx)*(max.size[sp]-7.5)/nx
  dbh.my.sp <- dbh.my.sp[1:nx]
  # size of the interval
  my.h <- dbh.my.sp[2] - dbh.my.sp[1]
  
  # arrange data, depending if originally "trees" is stored as a list or as a matrix
  if(type == "list"){
    index.trees <- GetIndex(trees)
    
    # which plots belong to that scenario
    my.trees <- index.trees[match(map$ID,index.trees[,2]),1]
    my.trees <- na.omit(my.trees)
    
    # add distributions to the average
    count <- 0 
    for(j in my.trees){
      if(!is.null(trees[[j]]$num[[sp]])){
        count <- count + 1
        dbh.distribution <- dbh.distribution + trees[[j]]$num[[sp]]
      }
    }
    
  }else{
    
    # I want a particular species
    my.trees <- trees[,sp,]
    
    if(length(plot.index != 0)){
      if(plot.index[1] != 0){
        my.trees <- my.trees[plot.index,]
      }
    }
    
    # how many plots have trees
    count <- length(which(rowSums(x=my.trees,dims=1)!=0))
    #     
    #     print(paste("sp:",sp,"- trees estimated:",quadTrapezCpp_1(q = colSums(my.trees[which(rowSums(my.trees,dims=1)!=0),],na.rm=T,dims=1),h = my.h,nx = nx),"- plots: ",count))
    
    # get the average distribution -- from the plots with presence, therefore the which
    # in other words, I don't want the average dbh distribution from all plots, because
    # in a majority of plots there might be no individuals
    
    # if there are plots with presence
    if(count>0){
      # if there are more than one plot
      if(count>1){
        dbh.distribution <- colMeans(my.trees[which(rowSums(my.trees,dims=1)!=0),],na.rm=T,dims=1)
      }else{ # only one plot
        dbh.distribution <- my.trees[which(rowSums(my.trees,dims=1)!=0),]
      }
    }#else{ # no plots
    #  dbh.distribution <- numeric(nx)
    #}
  }# if type list or matrix
  
  # evaluate integral: we want the mean number of trees in each size interval
  for(j in 1:length(dbh.my.sp)-1){
    dbh.distribution[j] <- interval.quad.trapez(q=dbh.distribution[c(j,j+1)],h=my.h)
  }
  
  # if there are no individuals, fill with 0s
  dbh.distribution[which(is.na(dbh.distribution))] <- 0
  
  results <- data.frame(frequency = dbh.distribution)
  results$scenario <- scenario
  results$year <- year
  results$dbh <- dbh.my.sp
  results$species <- sp
  
  if(write.to.file){
    write.table(results,paste(path,"DBH_distribution_sp",sp,"_",scenario,"_",year,"_",sim.number,".csv",sep=""))
  }
  
  results
  
}

#######################################
#######################################

SaveAbundance <- function(map,
                          trees,
                          saplings,
                          tesauro,
                          year,
                          scenario,
                          SIM_NUMBER,
                          plot.distribution = T,
                          plot.abundance = T,
                          plot.richness = T,
                          plot.ba = T,
                          spain.df,
                          path="./resultados/",
                          goparallel=F){
  # adults data frame
  adults.abundance <- data.frame(ID = map$ID,
                                 X_UTM = map$X_UTM, 
                                 Y_UTM = map$Y_UTM, 
                                 adult.trees = integer(nrow(map)),
                                 sp.basal.area = integer(nrow(map)))
  # saplings data frame
  saplings.abundance <- data.frame(ID = map$ID,
                                   X_UTM = map$X_UTM, 
                                   Y_UTM = map$Y_UTM, 
                                   saplings = integer(nrow(map)))
  # general data frame
  total.abundance <- data.frame(ID = map$ID,
                                X_UTM = map$X_UTM, 
                                Y_UTM = map$Y_UTM, 
                                num.trees = integer(nrow(map)),
                                plot.basal.area = integer(nrow(map)),
                                dominant.sp = integer(nrow(map)),
                                num.species = integer(nrow(map)))
  
  ########################
  ########################
  
  # auxiliary structures
  dominant.sp <- matrix(0,nrow=nrow(map),ncol=3)
  colnames(dominant.sp) <- c("basal.area","sp","num.species")
  
  num.trees <- numeric(nrow(map))
  
  for(i.sp in 1:nrow(tesauro)){
#     if(goparallel){
    # adults.abundance$adult.trees <- apply(X = trees[,i.sp,], MARGIN = 1, FUN = quadTrapezCpp_1,h[i.sp],nx) #unlist(mclapply(X=trees,FUN=helper.get.adults,i.sp))
    adults.abundance$adult.trees <- apply(X = trees[[i.sp]][], MARGIN = 1, FUN = quadTrapezCpp_1,h[i.sp],nx) #unlist(mclapply(X=trees,FUN=helper.get.adults,i.sp))
    adults.abundance$sp.basal.area <- ba[,i.sp] #unlist(mclapply(X=trees,FUN=helper.get.spba,i.sp))
    saplings.abundance$saplings <- saplings[,i.sp] * 10000/(pi*25) #unlist(mclapply(X=saplings,FUN=helper.get.saplings,i.sp)) * 10000/(pi*25)
#     }else{
#       adults.abundance$adult.trees <- unlist(lapply(X=trees,FUN=helper.get.adults,i.sp))
#       adults.abundance$sp.basal.area <- unlist(lapply(X=trees,FUN=helper.get.spba,i.sp))
#       saplings.abundance$saplings <- unlist(lapply(X=saplings,FUN=helper.get.saplings,i.sp)) * 10000/(pi*25)
#     }
    dominant.sp[,"sp"] <- ifelse(dominant.sp[,"sp"] > adults.abundance$sp.basal.area,dominant.sp[,"sp"],i.sp)
    dominant.sp[,"basal.area"] <- ifelse(dominant.sp[,"basal.area"] > adults.abundance$sp.basal.area,dominant.sp[,"basal.area"],adults.abundance$sp.basal.area)
    dominant.sp[,"num.species"] <- ifelse(adults.abundance$adult.trees>0,dominant.sp[,"num.species"]+1,dominant.sp[,"num.species"])
    num.trees <- num.trees + adults.abundance$adult.trees
    
    write.table(x=adults.abundance,paste(path,"trees_sp",i.sp,"_",year,"_",scenario,"_",SIM_NUMBER,".csv",sep=""),sep=";",dec=".",row.names=F,append=F)
    write.table(x=saplings.abundance,paste(path,"saplings_sp",i.sp,"_",year,"_",scenario,"_",SIM_NUMBER,".csv",sep=""),sep=";",dec=".",row.names=F,append=F)

    # plots
    if(sum(adults.abundance$adult.trees > 0) > 0){
      # distribution
      if(plot.distribution){
        plot.distribution(spain.df,
                          map,
                          adults.abundance,
                          i.sp,
                          year,
                          scenario,
                          SIM_NUMBER,
                          path=paste(path,"images",sep=""),
                          fill_ground = "#909090",
                          fill_points = "#505050",
                          x_axis = "UTM X",
                          y_axis = "UTM Y",
                          file.name = "dist_Sp",
                          plot.extension = "png")
      }
      # abundance
      if(plot.abundance){
        plot.abundance(spain.df=spain.df,
                       map=map,
                       trees.df=adults.abundance,
                       sp=i.sp,
                       year=year,
                       scenario=scenario,
                       SIM_NUMBER=SIM_NUMBER,
                       path=paste(path,"images",sep=""),
                       fill_ground = "#909090",
                       fill_points = "#505050",
                       x_axis = "UTM X",
                       y_axis = "UTM Y",
                       file.name = "abund_Sp",
                       plot.extension = "png")
      }
    }
  }
#   if(goparallel){
#     total.abundance$plot.basal.area <- unlist(mclapply(X=trees,FUN=helper.get.ba))
#   }else{
#     total.abundance$plot.basal.area <- unlist(lapply(X=trees,FUN=helper.get.ba))
#   }
  total.abundance$plot.basal.area <- apply(X = ba,MARGIN = 1,FUN = sum)#sum()
  total.abundance$num.trees <- num.trees
  total.abundance$dominant.sp <- dominant.sp[,"sp"]
  total.abundance$num.species <- dominant.sp[,"num.species"]
  
  write.table(x=total.abundance,file=paste(path,"total_abundance_",year,"_",scenario,"_",SIM_NUMBER,".csv",sep=""),sep=";",dec=".",row.names=F,append=F)
      
  # richness
  if(plot.richness){
    plot.richness(spain.df,
                  map,
                  total.abundance,
                  year,
                  scenario,
                  SIM_NUMBER,
                  path=paste(path,"images",sep=""),
                  fill_ground = "#909090",
                  fill_points = "#505050",
                  x_axis = "UTM X",
                  y_axis = "UTM Y",
                  file.name = "richness_",
                  plot.extension = "png")
  }
  # basal area
  if(plot.ba){
    plot.basal.area(spain.df,
                    map,
                    total.abundance,
                    year,
                    scenario,
                    SIM_NUMBER,
                    path=paste(path,"images",sep=""),
                    fill_ground = "#909090",
                    fill_points = "#505050",
                    x_axis = "UTM X",
                    y_axis = "UTM Y",
                    file.name = "basal_area_",
                    plot.extension = "png")
  }
}

#######################################
#######################################

SaveAbundance.list <- function(map,
                               trees,
                               saplings,
                               tesauro,
                               year,
                               scenario,
                               SIM_NUMBER,
                               plot.distribution = T,
                               plot.abundance = T,
                               plot.richness = T,
                               plot.ba = T,
                               spain.df,
                               path="./resultados/",
                               goparallel=F,
                               other.points=NULL){
  # adults data frame
  adults.abundance <- data.frame(ID = map$ID,
                                 X_UTM = map$X_UTM, 
                                 Y_UTM = map$Y_UTM, 
                                 adult.trees = integer(nrow(map)),
                                 sp.basal.area = integer(nrow(map)))
  # saplings data frame
  saplings.abundance <- data.frame(ID = map$ID,
                                   X_UTM = map$X_UTM, 
                                   Y_UTM = map$Y_UTM, 
                                   saplings = integer(nrow(map)))
  # general data frame
  total.abundance <- data.frame(ID = map$ID,
                                X_UTM = map$X_UTM, 
                                Y_UTM = map$Y_UTM, 
                                num.trees = integer(nrow(map)),
                                plot.basal.area = integer(nrow(map)),
                                dominant.sp = integer(nrow(map)),
                                num.species = integer(nrow(map)))
  
  ########################
  ########################
  
  # auxiliary structures
  dominant.sp <- matrix(nrow=nrow(map),ncol=3)
  colnames(dominant.sp) <- c("basal.area","sp","num.species")
  dominant.sp[,] <- 0
  
  num.trees <- numeric(nrow(map))
  
  for(i.sp in 1:nrow(tesauro)){

      adults.abundance$adult.trees <- unlist(lapply(X=trees,FUN=helper.get.adults,i.sp))
      adults.abundance$sp.basal.area <- ba[,i.sp] #unlist(mclapply(X=trees,FUN=helper.get.spba,i.sp))
      saplings.abundance$saplings <- saplings[,i.sp] * 10000/(pi*25)

    dominant.sp[,"sp"] <- ifelse(dominant.sp[,"sp"] > adults.abundance$sp.basal.area,dominant.sp[,"sp"],i.sp)
    dominant.sp[,"basal.area"] <- ifelse(dominant.sp[,"basal.area"] > adults.abundance$sp.basal.area,dominant.sp[,"basal.area"],adults.abundance$sp.basal.area)
    dominant.sp[,"num.species"] <- ifelse(adults.abundance$adult.trees>0,dominant.sp[,"num.species"]+1,dominant.sp[,"num.species"])
    num.trees <- num.trees + adults.abundance$adult.trees
    
    write.table(x=adults.abundance,paste(path,"trees_sp",i.sp,"_",year,"_",scenario,"_",SIM_NUMBER,".csv",sep=""),sep=";",dec=".",row.names=F,append=F)
    write.table(x=saplings.abundance,paste(path,"saplings_sp",i.sp,"_",year,"_",scenario,"_",SIM_NUMBER,".csv",sep=""),sep=";",dec=".",row.names=F,append=F)
    
    # plots
    if(sum(adults.abundance$adult.trees > 0) > 0){
      # distribution
      if(plot.distribution){
        plot.distribution(spain.df,
                          map,
                          adults.abundance,
                          i.sp,
                          year,
                          scenario,
                          SIM_NUMBER,
                          path="./resultados/images",
                          fill_ground = "#909090",
                          fill_points = "#505050",
                          x_axis = "UTM X",
                          y_axis = "UTM Y",
                          file.name = "dist_Sp",
                          plot.extension = "png")
      }
      # abundance
      if(plot.abundance){
        plot.abundance(spain.df=spain.df,
                       map=map,
                       trees.df=adults.abundance,
                       sp=i.sp,
                       year=year,
                       scenario=scenario,
                       SIM_NUMBER=SIM_NUMBER,
                       path="./resultados/images",
                       fill_ground = "#909090",
                       fill_points = "#505050",
                       x_axis = "UTM X",
                       y_axis = "UTM Y",
                       file.name = "abund_Sp",
                       plot.extension = "png",
                       other.points = other.points)
      }
    }
  }
  total.abundance$plot.basal.area <- apply(X = ba,MARGIN = 1,FUN = sum)#sum()
  total.abundance$num.trees <- num.trees
  total.abundance$dominant.sp <- dominant.sp[,"sp"]
  total.abundance$num.species <- dominant.sp[,"num.species"]
  
  write.table(x=total.abundance,file=paste(path,"total_abundance_",year,"_",scenario,"_",SIM_NUMBER,".csv",sep=""),sep=";",dec=".",row.names=F,append=F)
  
}

#######################################
#######################################

plot.distribution <- function(spain.df,
                  map,
                  trees.df,
                  sp,
                  year,
                  scenario,
                  SIM_NUMBER,
                  path="./resultados/images",
                  fill_ground = "#909090",
                  fill_points = "#505050",
                  x_axis = "UTM X",
                  y_axis = "UTM Y",
                  file.name = "dist_Sp",
                  plot.extension = "eps"){
  
  # check spatial dataframe consistence
  spain.df <- spain.df[spain.df$piece == 1,]
  
  points <- data.frame(ID = map$ID,
                       X_UTM = map$X_UTM,
                       Y_UTM = map$Y_UTM,
                       presence = integer(nrow(map)))
  
  if("adult.trees" %in% names(trees.df)){
    points$presence <- ifelse(trees.df$adult.trees > 0,1,0)
    points <- subset(points,presence > 0)
  }else{
    points$presence <- ifelse(trees.df$saplings > 0,1,0)
    points <- subset(points,presence > 0)
  }
  
  
  abc <- ggplot(spain.df)  
  abc <- abc + aes(long,lat) + xlab(x_axis) + ylab(y_axis)
  abc <- abc + geom_path(color="white") + coord_equal()
  abc <- abc + geom_polygon(fill=fill_ground)
  abc <- abc + geom_point(data=points,aes(x=X_UTM,y=Y_UTM),colour=fill_points,size=0.5)
  abc <- abc + theme(legend.position = "none") 
  
  ggsave(plot=abc,
         filename=paste(path,"/",file.name,"_",sp,"_",scenario,"_",year,"_v",SIM_NUMBER,".",plot.extension,sep=""),
         width=9,
         height=7,
         dpi=300)
  
}

#######################################
#######################################

plot.abundance <- function(spain.df,
                           map,
                           trees.df,
                           sp,
                           year,
                           scenario,
                           SIM_NUMBER,
                           path="./resultados/images",
                           fill_ground = "#909090",
                           fill_points = "#505050",
                           x_axis = "UTM X",
                           y_axis = "UTM Y",
                           file.name = "abund_Sp",
                           plot.extension = "eps",
                           other.points = NULL){
  # check spatial dataframe consistence
  spain.df <- spain.df[spain.df$piece == 1,]
  
  points <- data.frame(ID = map$ID,
                       X_UTM = map$X_UTM,
                       Y_UTM = map$Y_UTM,
                       abundance = integer(nrow(map)))
  
  if("adult.trees" %in% names(trees.df)){
    points$abundance <- trees.df$adult.trees
    points <- subset(points,abundance > 1)
  }else{
    points$abundance <- trees.df$saplings
    points <- subset(points,abundance > 0)
  }
  
  fun.cat <- function(x){
    if(x <= 50) cat <- 1
    if(x > 50 & x <= 100) cat <- 2
    if(x > 100) cat <- 3
    return(cat)
  }
  
  points$category <- sapply(points$abundance,fun.cat)
  points$category <- as.factor(points$category)
  
  abc <- ggplot(spain.df)  
  abc <- abc + aes(long,lat) + xlab(x_axis) + ylab(y_axis)
  abc <- abc + geom_path(color="black") + coord_equal()
  abc <- abc + geom_polygon(fill="white")
  # plot abundance with a color gradient
  abc <- abc + geom_point(data=points,aes(x=X_UTM,y=Y_UTM,color=category),size=0.00001) + scale_color_manual(values=c("#CCCCCC", "#999999", "#666666"),name="Adults abundance",labels=c("< 50 ind/ha", "50 - 100 ind/ha", "> 100 ind/ha")) #scale_color_gradient(low="grey",high="black",name = "Abundance")
  # do we have specific points to plot?
  
  if(!(is.null(other.points))){
    abc <- abc + geom_point(data=other.points,aes(x=X_UTM,y=Y_UTM,shape=category),size=3) + scale_shape_manual(values=c(16,17),name="Centroids",labels=c("adults","saplings")) 
  }# if points
  
  abc <- abc + guides(colour = guide_legend(override.aes = list(size=3)))
  abc <- abc + theme()#theme(legend.position = "none") 
  
  ggsave(plot=abc,filename=paste(path,"/",file.name,"_",sp,"_",scenario,"_",year,"_",SIM_NUMBER,".",plot.extension,sep=""),
         width=9,
         height=7,
         dpi=300)
}

#######################################
#######################################

  plot.richness <- function(spain.df,
                map,
                 summary.df,
                 year,
                 scenario,
                 SIM_NUMBER,
                 path="./resultados/images",
                 fill_ground = "#909090",
                 fill_points = "#505050",
                 x_axis = "UTM X",
                 y_axis = "UTM Y",
                 file.name = "richness",
                 plot.extension = "eps"){
    
    # check spatial dataframe consistence
    spain.df <- spain.df[spain.df$piece == 1,]
    
    points <- data.frame(ID = map$ID,
                         X_UTM = map$X_UTM,
                         Y_UTM = map$Y_UTM,
                         richness = integer(nrow(map)))
    
    points$richness <- summary.df$num.species
    points <- subset(points,richness > 0)
    levels <- length(unique(points$richness))
   
    abc <- ggplot(spain.df)  
    abc <- abc + aes(long,lat) + xlab(x_axis) + ylab(y_axis)
    abc <- abc + geom_path(color="black") + coord_equal()
    abc <- abc + geom_polygon(fill="white")
    # plot richness with a color gradient
    abc <- abc + geom_point(data=points,aes(x=X_UTM,y=Y_UTM,color=factor(richness)),size=0.5) + scale_color_brewer(palette = "Greens",name="Species Richness")#scale_color_manual(values = my.palette)
    
    ggsave(plot=abc,filename=paste(path,"/",file.name,"_",scenario,"_",year,"_v",SIM_NUMBER,".",plot.extension,sep=""),
           width=9,
           height=7,
           dpi=300)
    
}

#######################################
#######################################

plot.basal.area <- function(spain.df,
                map,
                  summary.df,
                  year,
                  scenario,
                  SIM_NUMBER,
                  path="./resultados/images",
                  fill_ground = "#909090",
                  fill_points = "#505050",
                  x_axis = "UTM X",
                  y_axis = "UTM Y",
                  file.name = "ba",
                  plot.extension = "eps"){
  # check spatial dataframe consistence
  spain.df <- spain.df[spain.df$piece == 1,]
  
  points <- data.frame(ID = map$ID,
                       X_UTM = map$X_UTM,
                       Y_UTM = map$Y_UTM,
                       richness = integer(nrow(map)))
  
  points$basal.area <- summary.df$plot.basal.area
  points <- subset(points,basal.area > 0)
  
  abc <- ggplot(spain.df)  
  abc <- abc + aes(long,lat) + xlab(x_axis) + ylab(y_axis)
  abc <- abc + geom_path(color="black") + coord_equal()
  abc <- abc + geom_polygon(fill="white")
  # plot basal arean with a color gradient
  abc <- abc + geom_point(data=points,aes(x=X_UTM,y=Y_UTM,color=basal.area),size=0.5) + scale_color_gradient(low="grey",high="black",name = "Basal Area")
  
  
  ggsave(plot=abc,filename=paste(path,"/",file.name,"_",scenario,"_",year,"_v",SIM_NUMBER,".",plot.extension,sep=""),
         width=9,
         height=7,
         dpi=300)
}

#######################################
#######################################

plot.gains.losses <- function(spain.df,
                  map,
                  trees.df.beginning,
                  trees.df.end,
                  sp,
                  year,
                  scenario,
                  SIM_NUMBER,
                  path="./resultados/images",
                  fill_ground = "#909090",
                  fill_points = "#505050",
                  x_axis = "UTM X",
                  y_axis = "UTM Y",
                  file.name = "gains_losses_Sp",
                  plot.extension = "eps"){
  
  # check spatial dataframe consistence
  spain.df <- spain.df[spain.df$piece == 1,]
  
  points.gained <- data.frame(ID = map$ID,
                       X_UTM = map$X_UTM,
                       Y_UTM = map$Y_UTM,
                       gained = logical(nrow(map)))
  
  points.lost <- data.frame(ID = map$ID,
                              X_UTM = map$X_UTM,
                              Y_UTM = map$Y_UTM,
                            lost = logical(nrow(map)))
  
  id.beginning <- unique(trees.df.beginning$ID)
  id.end <- unique(trees.df.end$ID)
  id.gained <- id.end[(is.na(match(id.end,id.beginning)))]
  id.lost <- id.beginning[is.na(match(id.beginning,id.end))]
  
  points.gained$gained <- sapply(points.gained$ID,function(x) return(x %in% id.gained))
  points.lost$lost <- sapply(points.lost$ID,function(x) return(x %in% id.lost))
  
  abc <- ggplot(spain.df)  
  abc <- abc + aes(long,lat) + xlab(x_axis) + ylab(y_axis)
  abc <- abc + geom_path(color="white") + coord_equal()
  abc <- abc + geom_polygon(fill=fill_ground)
  # plot gains with a color 
  abc <- abc + geom_point(data=points.gained,aes(x=X_UTM,y=Y_UTM),colour="green")
  # plot losses with another color 
  abc <- abc + geom_point(data=points.lost,aes(x=X_UTM,y=Y_UTM),colour="red")
  abc <- abc + theme(legend.position = "none") 
  
  ggsave(plot=abc,filename=paste(path,"/",file.name,"_",sp,"_",scenario,"_",year,"_v",SIM_NUMBER,".",plot.extension,sep=""))
}

#######################################
#######################################

# plot a temporal series of type "abundance", "basal area", or "both" 

plotTimeSeries <- function(type = "both",
                           year.list,
                           sp,
                           scenario,
                           SIM_NUMBER,
                           plot.ids = NULL,
                           path = "/home/david/CREAF/articulo/IPM/resultados"){
  
  abundance.list <- data.frame(year = sort(rep(year.list,length(sp))), sp = rep(sp,length(year.list)), abundance = 0, basal.area = 0)
  
  # load data
  for(i.year in 1:length(year.list)){
    for(i.sp in 1:length(sp)){
      
        my.file <- read.table(file = paste("/home/david/CREAF/articulo/IPM/resultados/trees_sp",sp[i.sp],"_",year.list[i.year],"_",scenario,"_",SIM_NUMBER,".csv",sep=""),header = T,sep = ";",dec = ".")
        
        # which plots
        my.file <- subset(my.file,adult.trees > 0 & sp.basal.area > 0)
        if(!(is.null(plot.ids))){
          my.file <- subset(my.file,ID %in% plot.ids)
        }

        abundance.list$abundance[abundance.list$year == year.list[i.year] & abundance.list$sp == sp[i.sp]] <- mean(my.file$adult.trees)
        abundance.list$basal.area[abundance.list$year == year.list[i.year] & abundance.list$sp == sp[i.sp]] <- mean(my.file$sp.basal.area)
        
    }# for i.sp
  }# for i.year
  
  sp.names <- character(length(sp))
  
  for(i.sp in 1:length(sp)){
    if(sp[i.sp] == 1){
      sp.names[i.sp] <- "Conifer"
    }else if(sp[i.sp] == 2){
      sp.names[i.sp] <- "Deciduous"
    }else if(sp[i.sp] == 3){
      sp.names[i.sp] <- "Fagus sylvatica"
    }else if(sp[i.sp] == 4){
      sp.names[i.sp] <- "Juniperus thurifera"
    }else if(sp[i.sp] == 5){
      sp.names[i.sp] <- "Pinus halepensis"
    }else if(sp[i.sp] == 6){
      sp.names[i.sp] <- "Pinus nigra"
    }else if(sp[i.sp] == 7){
      sp.names[i.sp] <- "Pinus pinaster"
    }else if(sp[i.sp] == 8){
      sp.names[i.sp] <- "Pinus pinea"
    }else if(sp[i.sp] == 9){
      sp.names[i.sp] <- "Pinus sylvestris"
    }else if(sp[i.sp] == 10){
      sp.names[i.sp] <- "Pinus uncinata"
    }else if(sp[i.sp] == 11){
      sp.names[i.sp] <- "Quercus faginea"
    }else if(sp[i.sp] == 12){
      sp.names[i.sp] <- "Quercus ilex"
    }else if(sp[i.sp] == 13){
      sp.names[i.sp] <- "Quercus pyrenaica"
    }else if(sp[i.sp] == 14){
      sp.names[i.sp] <- "Quercus robur/petraea"
    }else if(sp[i.sp] == 15){
      sp.names[i.sp] <- "Quercus suber"
    }else if(sp[i.sp] == 16){
      sp.names[i.sp] <- "Sclerophyllous"
    }
  }# for i.sp
  
  abund.plot <- ggplot(abundance.list)
  abund.plot <- abund.plot + geom_line(aes(x = year, y = abundance, linetype = as.factor(sp)))
  abund.plot <- abund.plot + xlab("year") + ylab("number of trees")
  abund.plot <- abund.plot + theme_JDGC() + scale_linetype_discrete(name = "Species",
                                                                    labels = sp.names)
  ba.plot <- ggplot(abundance.list)
  ba.plot <- ba.plot + geom_line(aes(x = year, y = basal.area, linetype = as.factor(sp)))
  ba.plot <- ba.plot + xlab("year") + ylab("basal area")
  ba.plot <- ba.plot + theme_JDGC() + scale_linetype_discrete(name = "Species",
                                                                    labels = sp.names)
  
  if(type == "both"){
    
    tt <- capture.output(cat(sp,sep="_"))
    ggsave(filename = paste(path,"/images/abund_sp_",tt,"_",year.list[1],"_",year.list[length(year.list)],"_",scenario,"_",SIM_NUMBER,".png",sep=""),
           plot = abund.plot,
           width = 7,
           height = 5,
           dpi = 300)
    
    ggsave(filename = paste(path,"/images/basal_area_sp_",tt,"_",year.list[1],"_",year.list[length(year.list)],"_",scenario,"_",SIM_NUMBER,".png",sep=""),
           plot = ba.plot,
           width = 7,
           height = 5,
           dpi = 300)
    
    plot.list <- list()
    plot.list[[1]] <- abund.plot
    plot.list[[2]] <- ba.plot
    plot.list
    
  }else if(type == "abundance"){
    tt <- capture.output(cat(sp,sep="_"))
    ggsave(filename = paste(path,"/images/abund_sp_",tt,"_",year.list[1],"_",year.list[length(year.list)],"_",scenario,"_",SIM_NUMBER,".png",sep=""),
           plot = abund.plot,
           width = 7,
           height = 5,
           dpi = 300)
    abund.plot
  }else{
    tt <- capture.output(cat(sp,sep="_"))
    ggsave(filename = paste(path,"/images/basal_area_sp_",tt,"_",year.list[1],"_",year.list[length(year.list)],"_",scenario,"_",SIM_NUMBER,".png",sep=""),
           plot = ba.plot,
           width = 7,
           height = 5,
           dpi = 300)
    ba.plot
  }
  
}

#######################################

GetSuitablePlots <- function(map,
                             temp.suitable,
                             ba,
                             #adult.trees,
                             #trees.index,
                             distancias,
                             sp,
                             max.dist,
                             verbose=FALSE,
                             goparallel=F){
  
  suitable <- temp.suitable
  
  #presence of the species
  #sum.ntrees <- apply(X = adult.trees[[sp]][],MARGIN = 1,FUN = sum) #apply(X = trees[,sp,],MARGIN = 1,FUN = sum)
  #presence.mysp <- ifelse(sum.ntrees>0,T,F)
  presence.mysp <- ifelse(ba[,sp]>0,T,F)
  
#   if(goparallel){
#     presence.mysp <- ifelse(trees[,sp,])
#     presence.mysp <- unlist(mclapply(trees,helper.get.presence,sp))
#   }else{
#     presence.mysp <- unlist(lapply(trees,helper.get.presence,sp))
#   }
  presence.ids <- map$ID[which(presence.mysp)]
  suitable <- suitable[which(!presence.mysp),]
  
  # which plots have lower basal area than threshold
  suitable <- suitable[suitable$sp_age==0,]
  suitable <- suitable[which(suitable$plot_basal_area < BA_threshold$perc_95[sp]),]
  
  #which have nearby presence of the species
  index <- match(suitable$ID,distancias[1,])
  
  if(length(index)>1){
    
    if(max.dist < 1414){
      close <- lapply(index,function(x) c(distancias[,x]$d1000))
    }else if(max.dist < 2000){
      close <- lapply(index,function(x) c(distancias[,x]$d1000, distancias[,x]$d1414))
    }else if(max.dist < 2236){
      close <- lapply(index,function(x) c(distancias[,x]$d1000, distancias[,x]$d1414, distancias[,x]$d2000))
    }else{
      close <- lapply(index,function(x) c(distancias[,x]$d1000, distancias[,x]$d1414, distancias[,x]$d2000, distancias[,x]$d2236))
    }
    
    if(class(close) == "list"){
      suitable$suitable <- sapply(close,function(x) sum(match(x,presence.ids,nomatch=0))>0)
      
      index <- which(suitable$suitable == TRUE)
      close <- close[index]
      close <- relist(match(unlist(close),map$ID),skeleton=close)
#       close <- relist(match(unlist(close),trees.index[,2]),skeleton=close)
      
      # clean up the index of nearby plots: we may find NAs, since we removed some rows from the map (invalid land use, plantation/riverside species...)
      close <- sapply(close,function(x) as.integer(na.omit(x)))
      
      suitable <- suitable[index,]
      if(verbose) print(paste(date()," - sp ",sp," - analyzing ",nrow(suitable)," plots...",sep=""))
      
      if(nrow(suitable)>0){
        # it's annoying not being able to do sum(trees[[close[i]]]$ba[[sp]]) sum(sapply(...)) is the workaround
        
        for (i in 1:nrow(suitable)) suitable$neigh_ba[i] <- sum(sapply(X=unlist(close[i]),FUN=function(x) ba[x,sp]))#trees[[x]]$ba[[sp]]))
        suitable <- suitable[,c("ID","plot_basal_area","neigh_ba")]
      }
    } else {
      suitable <- NULL
    }
  } else {
    suitable <- NULL
  }# if > 1
  suitable <- subset(suitable,neigh_ba > 20)
  return (suitable)
}

#######################################
#######################################

GetSuitablePlots.list <- function(map,
                                  temp.suitable,
                                  trees,
                                  ba,
                                  trees.index,
                                  distancias,
                                  sp,
                                  verbose=FALSE,
                                  max.dist=1000,
                                  goparallel=F){
  
  suitable <- temp.suitable
  
  #presence of the species
  presence.mysp <- ba[,sp]>0
  presence.ids <- map$ID[which(presence.mysp)]
  suitable <- suitable[which(!presence.mysp),]
  
  # which plots have lower basal area than threshold
  suitable <- suitable[which(suitable$plot_basal_area < BA_threshold$perc_95[sp]),]
  
  #which have nearby presence of the species
  index <- match(suitable$ID,distancias[1,])
  
  if(length(index)>1){
    
    if(max.dist < 1414){
      close <- lapply(index,function(x) c(distancias[,x]$d1000))
    }else if(max.dist < 2000){
      close <- lapply(index,function(x) c(distancias[,x]$d1000, distancias[,x]$d1414))
    }else if(max.dist < 2236){
      close <- lapply(index,function(x) c(distancias[,x]$d1000, distancias[,x]$d1414, distancias[,x]$d2000))
    }else{
      close <- lapply(index,function(x) c(distancias[,x]$d1000, distancias[,x]$d1414, distancias[,x]$d2000, distancias[,x]$d2236))
    }
    
    if(class(close) == "list"){
      suitable$suitable <- sapply(close,function(x) sum(match(x,presence.ids,nomatch=0))>0)
      
      index <- which(suitable$suitable == TRUE)
      close <- close[index]
      close <- relist(match(unlist(close),trees.index[,2]),skeleton=close)
      
      # clean up the index of nearby plots: we may find NAs, since we removed some rows from the map (invalid land use, plantation/riverside species...)
      close <- sapply(close,function(x) as.integer(na.omit(x)))
      
      suitable <- suitable[index,]
      if(verbose) print(paste(date()," - sp ",sp," - analyzing ",nrow(suitable)," plots...",sep=""))
      
      if(nrow(suitable)>0){
        # it's annoying not being able to do sum(trees[[close[i]]]$ba[[sp]]) sum(sapply(...)) is the workaround
        for (i in 1:nrow(suitable)) suitable$neigh_ba[i] <- sum(sapply(X=unlist(close[i]),FUN=function(x) ba[x,sp]))
        suitable <- suitable[,c("ID","plot_basal_area","neigh_ba")]
      }
    } else {
      suitable <- NULL
    }
  } else {
    suitable <- NULL
  }# if > 1
  return (suitable)
}

#######################################
#######################################

UTMDistance2 <- function(stand1,X,Y){
# Calculates the euclidean distance between two stands with UTM coordinates
#
# Args: stand1, X,Y. X and Y are vectors of x and y coordinates 
# Returns: the distance between the two stands (same units as original coord.)
#
   return (round(sqrt((stand1$X_UTM - X)^2 + (stand1$Y_UTM - Y)^2)))
}

#######################################
#######################################

Neighbours <- function(idparcela,map,mindist,maxdist){
  # returns id and coordinates of all neighbours of a plot closer than neighdist
  #
  #
  coordX <- map$X_UTM
  coordY <- map$Y_UTM
  dist <- UTMDistance2(map[map$ID == idparcela,],coordX,coordY)
  neigh <- map[which(dist > mindist & dist < maxdist),]
  return(neigh)
}

#######################################
#######################################

BasalArea <- function(pmfile,idparcela,sp = 0){
  # Calculates the total basal area of a given stand.
  # If sp > 0, only the basal area of the selected sp.
  # if sp == 0, all species
  #
  # Args: pmfile - file with all adult trees of all stands
  #       idparcela - the id of the stand to be calculated, or vector of IDs
  #       IFN - which IFN do we want calculations for
  # Returns: integer with the total basal area (m2/ha)
  #
  
  # which sp? if sp == 0, all of them. otherwise, pick the selected
  if (sp != 0){    
    pmfile <- subset(pmfile,num.sp == sp)
  }
  
  if(length(idparcela) == 1){
    
    basal.area <- 0
    temp.trees <- subset(pmfile, ID == idparcela)
    
    # if there are adult trees
    if (nrow(temp.trees) > 0){
      # basal area of all trees
      #temp.trees$basal.area <- pi*temp.trees$factor*(temp.trees$dbh/200)^2 
      #temp.trees$basal.area <- pi/4*(temp.trees$DN_IFN2*0.01)^2*temp.trees$FACTOR_IFN2
      # sum
      basal.area <- sum(temp.trees$factor*(temp.trees$dbh/200)^2)*pi
      
    } else basal.area <- 0
    
  }else{ # if idparcela > 1
    
    # the commented approach does not fill the stands that have zero basal area
    
    # temp.trees <- match(pmfile$IDPARCELA,idparcela)
    # temp.trees <- pmfile[!is.na(temp.trees),]
    #     
    # if(IFN == 2){
    #   temp.trees$basal.area <- pi/4 * (temp.trees$DN_IFN2 * 0.01)^2 * temp.trees$FACTOR_IFN2
    # }else{
    #   temp.trees$basal.area <- pi/4 * (temp.trees$DN_IFN3 * 0.01)^2 * temp.trees$FACTOR_IFN3
    # }
    #     
    # # sum all the basal areas by IDPARCELA
    # # aggregate returns a data frame with custom column names
    # basal.area <- aggregate(list(BASAL.AREA = temp.trees$basal.area),by=list(IDPARCELA=temp.trees$IDPARCELA),sum)
    # # equivalent is using ddply (from package plyr):
    # # basal.area <- ddply(temp.trees,"IDPARCELA",summarise, BASAL.AREA = sum(basal.area))
    #     
    
    # this is far less efficient, but works with all stands
    
    basal.area <- numeric(length(idparcela))
    for(i in 1:length(basal.area)){
      
      temp.trees <- subset(pmfile,ID == idparcela[i])
      
      if (nrow(temp.trees) > 0){
        temp.trees$basal.area <- pi*temp.trees$factor*((temp.trees$dbh/200)^2)
        #temp.trees$basal.area <- pi/4*(temp.trees$DN_IFN3*0.01)^2*temp.trees$FACTOR_IFN3
        basal.area[i] <- sum(temp.trees$basal.area)
      }
      else basal.area[i] <- 0  
    }# for basal.area
  }# if-else idparcela
  return (basal.area)
}

#######################################
#######################################

SurroundingBasalArea <- function(pmfile,mapfile,idparcela,mindist,maxdist,sp){
  # Calculates the basal area of a given species in surrounding stands, within a distance interval
  # if sp == 0, all 21 sp. are included
  #  
  # Args: pmfile - file with all adult trees of all stands
  #       mapfile - file with coordinates, necessary for distance calculations
  #       idparcela - the id of the central stand, or vector of IDs
  #       sp - species of which basal area is to be calculated
  #       mindist - minimum distance (m)
  #       maxdist - maximum distance (m)
  #       IFN - which IFN do we want calculations for
  # Returns: numeric with the total basal area; -1 if error with idparcela or mapfile
  #  
  
  Xcoord <- mapfile$X_UTM
  Ycoord <- mapfile$Y_UTM
  
  if(length(idparcela) == 1){
    
    total <- 0
    
    # there is only one central stand
    if (length(which(mapfile$ID == idparcela)) == 1){
      
      my.stand.coord <- mapfile[mapfile$ID == idparcela,]
      #distances to all stands
      dist <- UTMDistance2(my.stand.coord,Xcoord,Ycoord)
      #between mindist and maxdist
      temp <- which(dist >= mindist & dist <= maxdist)
      
      if(length(temp) > 0){
        total <- sum(BasalArea(pmfile,mapfile$ID[temp],sp)) 
      } # else keep total = 0
    }else total <- -1 # if not one central stand, return -1
    
  }else{ # if idparcela > 1
    
    total <- numeric(length(idparcela))
    
    for(j in 1:length(total)){
      # there is only one central stand
      if (length(which(mapfile$ID == idparcela[j])) == 1){
        
        my.stand.coord <- mapfile[mapfile$ID == idparcela[j],]
        #distances to all stands
        dist <- UTMDistance2(my.stand.coord,Xcoord,Ycoord)
        #between mindist and maxdist
        temp <- which(dist >= mindist & dist <= maxdist)
        
        if(length(temp) > 0){
          # sum of all the basal areas of "temp" stands
          total[j] <- sum(BasalArea(pmfile,mapfile$ID[temp],sp)) 
        } # else keep total[j] = 0
      } else total[j] <- -1 # if not one central stand, return -1
    }# for length total
  }#if-else idparcela
  
  return (total)
}

#######################################
#######################################

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#######################################
#######################################
# auxiliary plotting functions

stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, colour="red", geom=geom, width=0.2, ...)
}

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N    = length2(xx[,col], na.rm=na.rm),
                      mean = mean   (xx[,col], na.rm=na.rm),
                      sd   = sd     (xx[,col], na.rm=na.rm)
                   )
                 },
                 measurevar,
                 na.rm
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean"=measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#######################################
#######################################

# this function (from http://stackoverflow.com/questions/13297155/add-floating-axis-labels-in-facet-wrap-plot?lq=1)
# adjusts the x-axis on facet_wrap with uneven number of columns

# pos - where to add new labels
# newpage, vp - see ?print.ggplot
facetAdjust <- function(x, pos = c("up", "down"), 
                        newpage = is.null(vp), vp = NULL)
{
  # part of print.ggplot
  ggplot2:::set_last_plot(x)
  if(newpage)
    grid.newpage()
  pos <- match.arg(pos)
  p <- ggplot_build(x)
  gtable <- ggplot_gtable(p)
  # finding dimensions
  dims <- apply(p$panel$layout[2:3], 2, max)
  nrow <- dims[1]
  ncol <- dims[2]
  # number of panels in the plot
  panels <- sum(grepl("panel", names(gtable$grobs)))
  space <- ncol * nrow
  # missing panels
  n <- space - panels
  # checking whether modifications are needed
  if(panels != space){
    # indices of panels to fix
    idx <- (space - ncol - n + 1):(space - ncol)
    # copying x-axis of the last existing panel to the chosen panels 
    # in the row above
    gtable$grobs[paste0("axis_b",idx)] <- list(gtable$grobs[[paste0("axis_b",panels)]])
    if(pos == "down"){
      # if pos == down then shifting labels down to the same level as 
      # the x-axis of last panel
      rows <- grep(paste0("axis_b\\-[", idx[1], "-", idx[n], "]"), 
                   gtable$layout$name)
      lastAxis <- grep(paste0("axis_b\\-", panels), gtable$layout$name)
      gtable$layout[rows, c("t","b")] <- gtable$layout[lastAxis, c("t")]
    }
  }
  # again part of print.ggplot, plotting adjusted version
  if(is.null(vp)){
    grid.draw(gtable)
  }
  else{
    if (is.character(vp)) 
      seekViewport(vp)
    else pushViewport(vp)
    grid.draw(gtable)
    upViewport()
  }
  invisible(p)
}

list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  
  all <- list.files(path, pattern, all.dirs,
                    full.names, recursive=FALSE, ignore.case)
  all[file.info(all)$isdir]
}

base_breaks_x <- function(x){
  b <- pretty(x)
  d <- data.frame(y=-Inf, yend=-Inf, x=min(b), xend=max(b))
  list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend)),
       scale_x_continuous(breaks=b))
}
base_breaks_y <- function(x){
  b <- pretty(x)
  d <- data.frame(x=-Inf, xend=-Inf, y=min(b), yend=max(b))
  list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend)),
       scale_y_continuous(breaks=b))
}

#######################################
#######################################

theme_JDGC <- function(base_size = 10, base_family = "Helvetica") {
  theme(
    line =               element_line(colour = "black", size = 0.5, linetype = 1,
                                      lineend = "butt"),
    rect =               element_rect(fill = "black", colour = "black", size = 0.5, linetype = 1),
    text =               element_text(family = base_family, face = "plain",
                                      colour = "black", size = base_size,
                                      hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9),
    axis.text =          element_text(size = rel(0.8), colour = "black"),
    strip.text =         element_text(size = rel(0.8), colour = "black"),
    
    axis.line =          element_line(colour = "black", size = 0.5, linetype = 1,
                                      lineend = "butt"),#element_blank(),
    axis.text.x =        element_text(vjust = 1,),
    axis.text.y =        element_text(hjust = 1),
    axis.ticks =         element_line(colour = "black", size = 0.2),
    axis.title =         element_text(colour = "black"),
    axis.title.x =       element_text(vjust = 0.5),
    axis.title.y =       element_text(angle = 90),
    axis.ticks.length =  unit(0.2, "lines"),
    axis.ticks.margin =  unit(0.2, "lines"),
    
    legend.background =  element_rect(colour = NA,fill="white"),
    legend.margin =      unit(0.2, "cm"),
    legend.key =         element_rect(fill = "white", colour = "white"),
    legend.key.size =    unit(1.5, "lines"),
    legend.key.height =  NULL,
    legend.key.width =   NULL,
    legend.text =        element_text(size = rel(0.9), colour = "black"),
    legend.text.align =  NULL,
    legend.title =       element_text(size = rel(1), face = "bold", hjust = 0, colour = "black"),
    legend.title.align = NULL,
    legend.position =    "right",
    legend.direction =   "vertical",
    legend.justification = "bottom",
    legend.box =         NULL,
    
    panel.background =   element_rect(fill = "white", colour = NA),
    panel.border =       element_rect(fill = NA, colour = "black"),
    panel.grid.major =   element_blank(),#element_line(colour = "grey20", size = 0.2),
    panel.grid.minor =   element_blank(),#element_line(colour = "grey5", size = 0.5),
    panel.margin =       unit(0.25, "lines"),
    
    strip.background =   element_rect(fill = "lightgrey", colour = "black"),
    strip.text.x =       element_text(size=base_size),
    strip.text.y =       element_text(angle = -90,size=base_size),
    
    plot.background =    element_rect(colour = "white", fill = "white"),
    plot.title =         element_text(size = rel(1.2)),
    plot.margin =        unit(c(1, 1, 0.5, 0.5), "lines"),
    
    complete = TRUE
  )
}

#######################################
#######################################

# look for an expression in a string

getExpression <- function(astring,exp.vector){
  exp <- 0
  for(i.exp in 1:length(exp.vector)){
    if(any(grepl(pattern = exp.vector[i.exp],x = astring))){
      exp <- exp.vector[i.exp]
    }
  }# for
  return(exp)
}

#######################################
#######################################

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by = "Size", decreasing=TRUE, head = TRUE, n = 10) {
  # based on postings by Petr Pikal and David Hinds to the r-help list in 2004
  # modified by: Dirk Eddelbuettel (http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session) 
  # I then gave it a few tweaks (show size as megabytes and use defaults that I like)
  # a data frame of the objects and their associated storage needs.
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size) / 10^6 # megabytes
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

#######################################
#######################################

# standardize vector
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# standard error
se <- function(x) sqrt(var(x)/length(x))

