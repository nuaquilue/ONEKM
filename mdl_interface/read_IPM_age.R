read.IPM.age <- function(ba.file){
	load("./mdl_interface/cell.id.interface.rdata")
	load(ba.file)
	load("./mdl_interface/output/land.rdata")
	IPM.forest.age <- matrix(0,nrow = nrow(ba),ncol = ncol(ba))

	for (i in 1:nrow(ba)){
		Medfire.id <- cell.id.interface[which(cell.id.interface[,"IPM.index"]==i), "Medfire.id"]
		for (j in 1:ncol(ba)){
			if (ba[i,j]==0){
				IPM.forest.age[i,j] <- 0				
			}
			else{
				IPM.forest.age[i,j] <- land$age[Medfire.id]
			}
		}
	}
 	save(IPM.forest.age, file="./mdl_interface/IPM.forest.age.rdata")
 }