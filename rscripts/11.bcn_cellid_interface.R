library(raster)
IFN.map <- read.table("ONEKM/IPM/MAP_BCN_v1.csv",header=T,sep=";",dec=".")


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
Medfire.raster <- raster("./ONEKM/Medfire/inputlyrs/asc/LCFspp_1km_31N-ETRS89.asc")
IPM.raster <- raster()
extent(IPM.raster)<-extent(Medfire.raster)
res(IPM.raster)<-c(1000,1000) #1km^2
crs(IPM.raster) <- CRS("+init=epsg:25831")
IPM.raster <-rasterize(IFN_coord_31UTM, IPM.raster, IFN.map$ID)

cell.id.interface.bcn <- data.frame( Medfire.id=1:ncell(Medfire.raster), IFN.id=IPM.raster[] )
cell.id.interface.bcn <- cell.id.interface.bcn[!is.na(cell.id.interface.bcn$Medfire.id),]
cell.id.interface.bcn <- cell.id.interface.bcn[!is.na(cell.id.interface.bcn$IFN.id),]