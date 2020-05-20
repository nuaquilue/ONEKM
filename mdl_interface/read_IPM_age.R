read.IPM.age <- function(study.area, orig.plots.age.file){
	load("./mdl/inputlyrs/rdata/land.rdata")
	load((paste0(work.path,"/mdl_interface/IPM_",study.area, "_map.rdata"))
	IPM.forest.age <- matrix(0,nrow = nrow(ba),ncol = ncol(ba))

	for (i in 1:nrow(ba)){
		Medfire.id <- map$Medfire.id[i]
		for (j in 1:ncol(ba)){
			if (ba[i,j]==0){
				IPM.forest.age[i,j] <- 0				
			}
			else{
				IPM.forest.age[i,j] <- land$age[Medfire.id]
			}
		}
	}
 	save(IPM.forest.age, file=orig.plots.age.file)
 }