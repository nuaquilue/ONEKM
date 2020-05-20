#####################################################################################################
## INPUT
## work.path: local directory that should contain folders: inputlyrs and IPM 
######################################################################################################
## OUTPUT
## Builds 1 .rdata:
## cell.id.interface: associates cell.id between IPM plots and Medfire grid cells
######################################################################################################
## NOTES
## This function assumes that Medfire projections are in  31N-ETRS89 and IPM are in 30N-ED50
######################################################################################################

make.cell.ID.interface <- function(work.path){

  ##IFN.map: IFN spain map
  setwd(work.path)
  IFN.map <- read.table(paste0(work.path,"/IPM/MAP_v6.csv"),header=T,sep=";",dec=".")

  
  ##Change X_UTM of IPM so they are in 31 N projection (Medfire) instead of 30 N 
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
  IPM.raster <-rasterize(IFN_coord_31UTM, IPM.raster, IFN.map$ID)

  ##create cell.id.interface dataframe:
  ## Medfire.id | IFN.id 
  cell.id.interface <- data.frame( Medfire.id=1:ncell(Medfire.raster), IFN.id=IPM.raster[] )
  cell.id.interface <-cell.id.interface[!is.na(Medfire.raster[]),]
  save(cell.id.interface, file=paste0(work.path,"/mdl_interface/cell.id.interface.rdata"))
 
}