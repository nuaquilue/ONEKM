select.valid.plots <- function(study.area){

	load("./mdl_interface/cell.id.interface.rdata")
	load("inputlyrs/rdata/land.rdata")
	if (study.area=="BCN"){
		map <- read.table(file="./IPM/MAP_BCN_v1.csv",header=T,sep=";",dec=".")
	} else if(study.area=="CAT"){
		map <- read.table(file="./IPM/MAP_CAT_v2.csv",header=T,sep=";",dec=".")
	} else print("to work with a study area different of CAT or BCN, it is needed to add csv file and modify function\n")
	
	map <- left_join(map, cell.id.interface, by=c("ID"="IFN.id")) %>% left_join(land[,c("cell.id","spp")], by=c("Medfire.id"="cell.id"))
	#elim_plots <- map[is.na(map$Medfire.id),] ##do we want to save them?
	map <- map[!is.na(map$Medfire.id),] ##Eliminates plots that are outside Medfire study area
	map <- map[map$spp<15,]
	save(map, file=paste0(work.path,"/mdl_interface/IPM_",study.area, "_map.rdata"))
  
}