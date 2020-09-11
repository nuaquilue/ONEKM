library(raster)
work.path <- "C:/Users/uriso/Desktop/Medfire_analysis"
####################
##ONEKM
####################
setwd(work.path)
fire.zones <- raster("./FireRegimeZones_31N-ETRS89_1km.asc")
##fire.zones <- raster("./FireRegimeZones_31N-ETRS89_1km.asc")
zones_poly <- fire.zones[]
zones_poly[which(is.na(zones_poly))]<- 80
limit_cells<-c()
for(index in 291:(length(zones_poly)-290)){
  
  if (index%%290!=0 & (index-1)%%290!=0){
    sum_neigh <- sum(zones_poly[index-1],zones_poly[index+1],
                       zones_poly[index+290],zones_poly[index-290])
    i <- zones_poly[index]
    if (sum_neigh!= (i*4)){
        limit_cells<-c(limit_cells, index)
      }
  }
  
  index<- index+1
}
poli<- fire.zones
poli[]<- NA
poli[limit_cells]<-1


tb_onekm <- ("../tb_onekm")
clim.scn <- "rcp45"
tif_files_onekm <- list.files(path=paste0(tb_onekm, "/", clim.scn), pattern="*.tif", full.names=TRUE, recursive=FALSE)

num_tif <- length(tif_files_onekm)

tif1<- raster(tif_files_onekm[1])
mean_tb_onekm <- tif1[]
for (i in 2:num_tif){
  tif_i<- raster(tif_files_onekm[i])
  mean_tb_onekm<- mean_tb_onekm+tif_i[]
}
mean_tb_onekm<- mean_tb_onekm/num_tif

num_fire_zones <- length(unique(fire.zones[]))-1
##times burnt only accounts for spp <=17, therefore we have to apply mask
##LCF<-raster("./LCFspp_1km_31N-ETRS89.asc")
##fire.zones[LCF[]>17]<- NA
#fire.zones[is.na(tif1[])]<- NA
years_raster <- fire.zones

for(i in 1:num_fire_zones){
  zone_i <- which(fire.zones[]==i) ##ID 
  len_zone_i <- length(zone_i)
  tb_zone_i <- mean_tb_onekm[zone_i]
  years_to_burn <- len_zone_i/sum(tb_zone_i, na.rm=TRUE)*100
  years_raster[zone_i] <- years_to_burn
}
years_raster[years_raster>900]<- 1000
#years_raster[limit_cells]<-1500
cuts<-c(50,100,200,300,400,600,900,1200)
pal <- colorRampPalette(c("red","yellow", "green"))
plot(years_raster, breaks=cuts, col=c(pal(8), "black"))
scalebar(100000, xy=c(390000, 4725000), type='bar', divs=4, label=c(0, 50, 100), below="Km", cex=0.5)

###############################
##ONEHAC
###############################
fire.zones <- raster("./FireRegimeZones_31N-ETRS89_100m.asc")
##fire.zones <- raster("./FireRegimeZones_31N-ETRS89_1km.asc")

tb_onehac <- ("../tb_onehac")
clim.scn <- "Scn_rcp45_noFF_2/lyr"
tif_files_onehac <- list.files(path=paste0(tb_onehac, "/", clim.scn), pattern="*.tif", full.names=TRUE, recursive=FALSE)

num_tif <- length(tif_files_onehac)

tif1<- raster(tif_files_onehac[1])
mean_tb_onehac <- tif1[]
for (i in 2:num_tif){
  tif_i<- raster(tif_files_onehac[i])
  mean_tb_onehac<- mean_tb_onehac+tif_i[]
}
mean_tb_onehac<- mean_tb_onehac/num_tif

num_fire_zones <- length(unique(fire.zones[]))-1
##times burnt only accounts for spp <=17, therefore we have to apply mask
##LCF<-raster("./LCFspp_1km_31N-ETRS89.asc")
##fire.zones[LCF[]>17]<- NA
#fire.zones[is.na(tif1[])]<- NA
years_raster <- fire.zones

for(i in 1:num_fire_zones){
  zone_i <- which(fire.zones[]==i) ##ID 
  len_zone_i <- length(zone_i)
  tb_zone_i <- mean_tb_onehac[zone_i]
  years_to_burn <- len_zone_i/sum(tb_zone_i, na.rm=TRUE)*100
  years_raster[zone_i] <- years_to_burn
  cat(paste0("analysing zone: ", i, "\n"))
}
years_raster[years_raster>900]<- 1000

cuts<-c(50,100,200,300,400,600,900,1200)
pal <- colorRampPalette(c("red","yellow", "green"))
plot(years_raster, breaks=cuts, col=pal(8))
##plot legend
# plot(1,2, xlim=c(0, 34), pch= 15,col= pal(8)[1], cex=2)
# for(i in 2:8){
#   points((4*i-3),2, pch= 15,col= pal(8)[i], cex=2)
# }
cuts <- c("0-50", "100-200", "200-300", "300-400", "400-600", "600-900", ">900")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend = cuts, pch=15, pt.cex=3, cex=1.5, bty='n',
       col = pal(8))
#mtext("Legend", at=0.2, cex=2)

##plotting difference
years_raster_onekm <- years_raster
load("./years_raster_onehac.rdata")
years_raster_onehac <- aggregate(years_raster,fact=c(10,10))
years_raster_onehac[which(is.na(years_raster_onekm[]))]<- NA


tb_diff <- years_raster_onekm #years_raster#
tb_diff[]<- years_raster_onehac[]/years_raster_onehac_v0[]*100-100 #years_raster[]/years_raster_v1_10[]*100-100 #
tb_diff_zones <- rep(0,num_fire_zones)
for(i in 1:num_fire_zones){
  zone_i <- which(fire.zones[]==i) ##ID 
  tb_diff_zones[i]<- mean(tb_diff[zone_i] , na.rm=TRUE)
  cat(paste0("analysing zone: ", i, "\n"))
}
plot(tb_diff)
mean(tb_diff[], na.rm=T)

##Fire distributions
##1hac
fires_onehac <- read.table(paste0(tb_onehac, "/", clim.scn,"/fires.txt") ,header = T)
sum_fires_onehac<-data.frame(matrix(NA, nrow = num_tif , ncol = 4))
for(run in 1:num_tif){
  irun <- which(fires_onehac$run==run)
  sum_fires_onehac[run,]<-c(sum(fires_onehac$atarget[irun]), sum(fires_onehac$aburnt.highintens[irun]),
                          sum(fires_onehac$aburnt.lowintens[irun]), sum(fires_onehac$rem[irun]))
  
}
##1km
fires_onekm <- read.table(paste0(tb_onekm, "/", clim.scn,"/fires_1.txt") ,header = T)
sum_fires_onekm<-data.frame(matrix(NA, nrow = num_tif , ncol = 4))
for(run in 1:5){
  irun <- which(fires_onekm$run==run)
  sum_fires_onekm[run,]<-c(sum(fires_onekm$atarget[irun]), sum(fires_onekm$aburnt.highintens[irun]),
                            sum(fires_onekm$aburnt.lowintens[irun]), sum(fires_onekm$rem[irun]))
  
}
fires_onekm <- read.table(paste0(tb_onekm, "/", clim.scn,"/fires_2.txt") ,header = T)
for(run in 1:5){
  irun <- which(fires_onekm$run==run)
  sum_fires_onekm[run+5,]<-c(sum(fires_onekm$atarget[irun]), sum(fires_onekm$aburnt.highintens[irun]),
                           sum(fires_onekm$aburnt.lowintens[irun]), sum(fires_onekm$rem[irun]))
  
}

library(ggplot2)
area_burnt <- data.frame(scale=rep(c("onehac", "onekm"), each=4),
                  mean_area=rep(c("atarget", "ahigh_int", "alow_int", "not_burnt"),2),
                  hac=c(mean(sum_fires_onehac[,1]), mean(sum_fires_onehac[,2]), mean(sum_fires_onehac[,3]),mean(sum_fires_onehac[,4]),
                        mean(sum_fires_onekm[,1]), mean(sum_fires_onekm[,2]), mean(sum_fires_onekm[,3]),mean(sum_fires_onekm[,4])))
ggplot(data=area_burnt, aes(x= factor(mean_area, level = c("atarget", "ahigh_int", "alow_int", "not_burnt")), y=hac, fill=scale)) +
geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=hac), vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

