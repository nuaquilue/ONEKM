library('plot.matrix')
ha_matrix <- matrix(nrow=20, ncol = 20)
ul <- sample(c(1,2,3,4), 100, replace = T, prob =c(0.6,0.2,0.1,0.1))
ur <- sample(c(1,2,3,4), 100, replace = T, prob =c(0.1,0.,0.2,0.1))
dl <- sample(c(1,2,3,4), 100, replace = T, prob =c(0.1,0.6,0.2,0.1))
dr <- sample(c(1,2,3,4), 100, replace = T, prob =c(0.1,0.3,0.6,0))
ha_matrix[1:10,1:10]<- ul
ha_matrix[1:10,11:20]<- ur
ha_matrix[11:20,1:10]<- dl
ha_matrix[11:20,11:20]<- dr
colors <- c("#FFFF99", "#CCFF66", "#99FF33", "#66CC00", "#339900")

plot(ha_matrix, xlab="", ylab="" , breaks=c(0.00,1.5,2.5,3.5,4.5),
     col=colors[-3], main="")
ul_km <- mode(ul)
ur_km <- mode(ur)
dl_km <- mode(dl)
dr_km <- mode(dr)
km_matrix<- matrix(nrow=2,ncol=2,c(1,2,2,3))
plot(km_matrix, xlab="", ylab="" , breaks=c(0.00,1.5,2.5,3.5,4.5),
     col=colors[-3], main="")
##plot rasters
study.area <- "BCN"
load(paste0("ONEKM/mdl_interface/IPM_",study.area, "_map.rdata"))

bcn_shape <- rgdal::readOGR(dsn="./IPM_10y",layer="bcn_UTM_v1")
berga <- rgdal::readOGR(dsn="./berga",layer="veg293_2014_ETRS89")

points <- data.frame(lat = map$X_UTM, lon = map$Y_UTM)
points <- SpatialPoints(points)
plot(bcn_shape)
plot(points, add= TRUE , pch=16, cex=0.01)

bcn_shape_2 <- bcn_shape
extent(bcn_shape_2)<- c(900000, 920000, 4620000,4640000)
points_2 <- crop(points, c(900000, 920000, 4620000,4640000))
#plot(map$X_UTM,map$Y_UTM, pch=17, bg="transparent", add=TRUE)

####Berga Raster

berga.raster <- raster()
extent(berga.raster)<-extent(berga)
res(berga.raster)<-c(1000,1000) #1km^2
crs(berga.raster) <- CRS("+init=epsg:25831")
berga.raster[]<-1

load("./ONEKM/Medfire/inputlyrs/rdata/coordinates.rdata")
points <- SpatialPoints(coord[coord$cell.id %in% map$Medfire.id[map$spp<=14],c(2,3)])
plot(berga.raster, col="w")
plot(points, add= TRUE , pch=16, cex=0.2)



#LCFspp_1km_31N-ETRS89
color <- c("chartreuse3", "darkolivegreen1", "darkseagreen", "forestgreen",
                 "olivedrab4", "darkslategrey", "blue4", "gold", "saddlebrown",
                 "sienna2", "palegoldenrod", "red3", "purple3", "grey70")
LCF <- raster( "./ONEKM/Medfire/inputlyrs/asc/LCFspp_1km_31N-ETRS89.asc")


santcel_LCT <- crop(LCF, extent(santcel_shp)+c(-2000,2000,-2000,2000))
santcel_LCT <- mask(santcel_LCT, santcel_shp)
#berga_LCT[berga_LCT>14]<- 25

cuts<-c(seq(0.5,14.5,1), 30)
plot(santcel_LCT, breaks=cuts, col=c(color, "gray"))
plot(santcel_shp, add= TRUE)

plot(points, add= TRUE , pch=16, cex=0.2)


##
capitals <- rgdal::readOGR(dsn="./remote_output/CapitalsComarca",layer="CapitalsComarca")
cat_shp <- rgdal::readOGR(dsn="./cat_shp",layer="MUC_TM")
# terrassa_shp <- cat_shp[cat_shp@data$MUNICIPI=="Terrassa",]
# montseny_shp <- cat_shp[cat_shp@data$MUNICIPI=="Montseny",]
santcel_shp <- cat_shp[cat_shp@data$MUNICIPI=="Sant Celoni",]
plot(santcel_shp)
over(points, as(santcel_shp,"SpatialPolygons"))
plot(points[which(!is.na(over(points, as(santcel_shp,"SpatialPolygons")))),], add= TRUE , pch=16, cex=0.2)


plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)

species <- c("phalepensis", "pnigra", "ppinea", "psylvestris", "ppinaster", "puncinata",
             "aalba", "qilex", "qsuber", "qfaginea", "qhumilis", "fsylvatica", "other")
legend("topleft", legend = c(  species, "shrub", "urban/crops" ), pch=15, pt.cex=1, cex=0.7, bty='n',
       col = c(color,"gray"), ncol =4)