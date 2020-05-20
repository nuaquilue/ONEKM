#####################################################################################################
## INPUT
## work.path: local directory that should contain folders: inputlyrs and IPM 
## map.csv: csv map of the study area in IPM
######################################################################################################
## OUTPUT
## Builds 2 .rdata:
## 1. cell.id.interface: associates cell.id and coordenates between IPM plots and Medfire grid cells
## 2. IPM.deleted.plots: IPM plots that cannot be associate to any Medfire grid cell
######################################################################################################
## NOTES
## This function assumes that Medfire projections are in  31N-ETRS89 and IPM are in 30N-ED50
######################################################################################################

make.cell.ID.interface <- function(work.path, map.csv){

  ##IFN.map: filters IPM Spain map (MAP_v6.csv) to only valid plots within target.area
  target.area.map <- read.table(file=paste0(work.path,"/IPM/",map.csv),header=T,sep=";",dec=".")
  IFN.map <- read.table(paste0(work.path,"/IPM/MAP_v6.csv"),header=T,sep=";",dec=".")
  IFN.map <- subset(IFN.map, valid)
  IFN.map <- subset(IFN.map, ID %in% target.area.map$ID)
  
  ##Change X_UTM of IPM so they are in 31 N projection (Medfire) instead of 30 N 
  ##If  Medfire projection changes, this needs to be changed.
  #Change coordinates of IPM so they are in 31 N projection (Medfire) instead of 30 N 
  ##If  Medfire projection changes, this needs to be changed.
  IFN_coord<-IFN.map[,c('X_UTM', 'Y_UTM')]
  colnames(IFN_coord)<- c("lon","lat")
  coordinates(IFN_coord) <-c("lon","lat")
  proj4string(IFN_coord) <- CRS("+init=epsg:23030")
  new_coord <-spTransform(IFN_coord,CRS("+init=epsg:25831"))
  IFN_coord_31UTM<- coordinates(new_coord)

  ##Creates raster file from IFN.map dataframe with Medfire rasters extension
  ##In order to associate IFN plots to same Medfire coordinates
  Medfire.raster <- raster("./inputlyrs/asc/LCFspp_1km_31N-ETRS89.asc")
  IPM.raster <- raster()
  extent(IPM.raster)<-extent(Medfire.raster)
  res(IPM.raster)<-c(1000,1000) #1km^2
  crs(IPM.raster) <- CRS("+init=epsg:25831")
  IPM.raster <-rasterize(IFN_coord_31UTM, IPM.raster, 1:nrow(IFN.map))

  ##create cell.id.interface dataframe:
  ## | IPM.index | Medfire.id | X.Medfire | Y.Medfire | IFN.id | X.IFN | Y.IFN |
  ## sorted by IPM.index
  cell.id.interface <- data.frame(IPM.index=IPM.raster[], Medfire.id=1:ncell(Medfire.raster), X.Medfire=coordinates(IPM.raster)[,1], Y.Medfire=coordinates(IPM.raster)[,2])
  cell.id.interface <-cell.id.interface[!is.na(cell.id.interface$IPM.index),]
  cell.id.interface <- cell.id.interface[order(cell.id.interface$IPM.index),]
  cell.id.interface$X.IFN <- IFN_coord_31UTM[cell.id.interface$IPM.indexIPM.index,1]
  cell.id.interface$Y.IFN <- IFN_coord_31UTM[cell.id.interface$IPM.index,2]
  cell.id.interface$IFN.id <- IFN.map$ID[cell.id.interface$IPM.index]
  save(cell.id.interface, file=paste0(work.path,"/mdl_interface/cell.id.interface.rdata"))
  
  ##look for IFN plots that disappear in the raster form
  ##this happens when 2 (or more) plots are in the same 1km square grid 
  ##in the raster coordinates. Since rasterize only stores one of them.
  del.plots.IPM.index <- (1:nrow(IFN.map))[!(1:nrow(IFN.map) %in% cell.id.interface$IPM.index)]
  del.plots.IFN.id <- IFN.map[del.plots.IPM.index,"ID"]
  IPM.deleted.plots <- data.frame(IPM.index= del.plots.IPM.index, IFN.id= del.plots.IFN.id)
  save(IPM.deleted.plots, file=paste0(work.path,"/mdl_interface/IPM.deleted.plots.rdata"))
  
}