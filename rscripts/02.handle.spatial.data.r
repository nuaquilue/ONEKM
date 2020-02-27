rm(list=ls())
library(raster)
library(rgdal)
## read raster
DEM <- raster("c:/work/MEDMOD/DataCLIM/DataSp/DEM1k.CAT_31N-ETRS89.asc")
## access raster values
dem <- DEM[]
## plot it
plot(DEM)
## assign cartographic projection: UTM31N-ETRS89
crs(DEM) <- CRS("+init=epsg:25831")
## compute slope and aspect and write raster
SLOPE <- terrain(DEM, opt='slope', unit='degrees', neighbors=8)
ASPECT <- terrain(DEM, opt='aspect', unit='degrees', neighbors=8)
writeRaster(SLOPE, "c:/work/MEDMOD/DataCLIM/DataSp/SlopeDegree_CAT1K.asc", format="ascii", overwrite=T)
writeRaster(ASPECT, "D:/MEDMOD/DataCLIM/DataSp/AspectDegree_CAT1K.asc")

## Read a SHP file: Vector of Line coast of Cataloniain UTM31N-ETRS89
COASTcat <- readOGR("c:/work/MEDMOD/DataCLIM/DataSp/Costa_ETRS89.shp")
## Plot it
plot(COASTcat)
## Assign cartographic projection
crs(COASTcat) <- CRS("+init=epsg:25831")
## Vector of world line coast in Long/Lat WGS84 
COASTall <- readOGR("c:/work/MEDMOD/DataCLIM/DataSp/coastline/ne_10m_coastline.shp")
## Cut to a smaller extent the above vector
COASTpeni <- crop(COASTall, extent(c(-15,5,35,45)))
## To change cartographic projection of a vector
COASTutm <- spTransform(COASTpeni, CRS("+init=epsg:25831"))
COAST <- crop(COASTutm, extent(DEM))


## Transform data frame to a raster 
res=1000
XORRADA <- rasterFromXYZ(cbind(coordinates(DEM), DEM[]*1000), res=c(res,res), digits=5)
## Assign new values to a raster
XORRADA[] <- ifelse(XORRADA[]>10^6, 9,8)
plot(XORRADA) 


## Change cartographic projection of a raster
tmin.file <- paste0("ClimHistAEMETsp/TNMM_", mdl,"_HISTORICAL_r1i1p1_", d,".txt")
TMIN <- raster(tmin.file)
crs(TMIN) <- CRS("+init=epsg:4326")  ## LatLong WGS84
TMINutm <- projectRaster(TMIN, res=10000, crs= CRS("+init=epsg:25831"))  ## UTM31N-ETRS89


## Read raster at 1Km, change resolution to 1ha, and change bounding box to the default one
## disaggregate vs aggregate
## Default extent of raster maps of Catalonia  
extCat <- extent(c(250000, 540000, 4480000, 4760000))
TEMP <- raster(paste0("inputlyrs/asc/", clim.scn, "/", decade, "/mnan.asc"))
TEMP <- disaggregate(TEMP, fact=c(10,10))
TEMP <- extend(TEMP, extCat)

