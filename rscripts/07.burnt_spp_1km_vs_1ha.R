work.path <- "C:/Users/uriso/Desktop/Medfire_analysis"
tb_onekm <- ("../tb_onekm")
clim.scn <- "rcp45"
fires_spp <- read.table(paste0(tb_onekm, "/", clim.scn,"/FiresSpp.txt"), header=T)
burnt_spp_1km <- matrix(nrow= 10, ncol=17)
for (run in 1:10){
  for (sp in 1:17){
    burnt_spp_1km[run,sp]<- sum(fires_spp[which(fires_spp$run==run & fires_spp$spp==sp),5])*100
  }
  
}

burnt_biom_1km <- matrix(nrow= 10, ncol=17)
for (run in 1:10){
  for (sp in 1:17){
    burnt_biom_1km[run,sp]<- sum(fires_spp[which(fires_spp$run==run & fires_spp$spp==sp),6])*100
  }
  
}

##100m
tb_onehac <- ("../tb_onehac")
fires_spp_100m <- read.table(paste0(tb_onehac, "/", clim.scn,"/FiresSpp.txt"), header=T)
burnt_spp_100m <- matrix(nrow= 10, ncol=17)
for (run in 1:10){
  for (sp in 1:17){
    burnt_spp_100m[run,sp]<- sum(fires_spp_100m[which(fires_spp_100m$run==run & fires_spp_100m$spp==sp),5])
  }
  
}

burnt_biom_100m <- matrix(nrow= 10, ncol=17)
for (run in 1:10){
  for (sp in 1:17){
    burnt_biom_100m[run,sp]<- sum(fires_spp_100m[which(fires_spp_100m$run==run & fires_spp_100m$spp==sp),6])
  }
  
}
species <- c("phalepensis", "pnigra", "ppinea", "psylvestris", "ppinaster", "puncinata",
             "aalba", "qilex", "qsuber", "qfaginea", "qhumilis", "fsylvatica", "other",
             "shrubland", "grass", "arable land", "crops")
##plot it
data <- (matrix(nrow=4, ncol=17))
data[1,]<-colMeans(burnt_spp_100m)
data[2,]<- colMeans(burnt_spp_1km)
data[3,]<- apply(burnt_spp_100m,2,sd)
data[4,]<- apply(burnt_spp_1km,2,sd)
num_spp_barplot <- barplot(data[c(1,2),], main="burnt land cover distribution", names.arg=species,
        ylab="ha",cex.axis =0.8,ylim=c(0,max(data[c(1,2),])+max(data[c(1,2),])*0.25) ,col=c("darkblue","cyan"),beside=TRUE, las=2,
        legend = c("MEDFIRE-100m", "MEDFIRE-1km"), args.legend = list(bty = "n", x = "top")) #legend = c("onehac", "onekm"),
arrows(num_spp_barplot, data[c(1,2),]-data[c(3,4),], num_spp_barplot, data[c(1,2),]+data[c(3,4),], angle=90, code=3, length=0.03 )
data <- matrix(nrow=2, ncol=13)
data[1,]<-colMeans(burnt_biom_100m)[1:13]
data[2,]<- colMeans(burnt_biom_1km)[1:13]
barplot(data, main="burnt land cover distribution",
        xlab="biomass", col=c("darkblue","red"),
        beside=TRUE) #legend = c("onehac", "onekm"),

diff_sp <- colMeans(burnt_spp_1km)/colMeans(burnt_spp_100m)*100-100
barplot(diff_sp, col="gray", main=("burnt 1km/burnt 1hac*100-100"),
        names.arg=species, las=2 )

diff_biom <- colMeans(burnt_biom_1km)[1:13]/colMeans(burnt_biom_100m)[1:13]*100-100
barplot(diff_sp)