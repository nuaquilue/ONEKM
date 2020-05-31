read.IPM.age <- function(ba, map){
	load("./Medfire/inputlyrs/rdata/land.rdata")
	#load(paste0(work.path,"/mdl_interface/IPM_",study.area, "_map.rdata"))
	IPM.forest.age <- matrix(0,nrow = nrow(ba),ncol = ncol(ba))

	for (i in 1:nrow(ba)){
		Medfire.id <- map$Medfire.id[i]
		for (j in 1:ncol(ba)){
			if (ba[i,j]==0){
				IPM.forest.age[i,j] <- 0				
			} else if (map$spp[i] !=14) {
				IPM.forest.age[i,j] <- land$age[which(land$cell.id == Medfire.id)]
			} else{
				IPM.forest.age[i,j] <- 10
			}
		}
	}
 	save(IPM.forest.age, file=orig.plots.age.file)
 }