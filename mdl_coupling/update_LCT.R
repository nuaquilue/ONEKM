update.LCF <- function(){
	load("./mdl_coupling/IPM_Medfire_Interface.rdata")
	load(ba.file)
	load("./mdl_coupling/output/land.rdata")
	for (i in 1:nrow(ba)){
		dom_species <- 0
		dom_abundance <- 0
		for (j in 1:length(newspecies)){
			if (ba[i,j]>dom_abundance){
				dom_abundance <- ba[i,j]
				dom_species <- j
			}
		}
		IPM_Medfire_Interface[i,"IPM_LCT"] <- dom_species
		if (dom_species!=0){
			IPM_Medfire_Interface[i,"IPM_LCT_Medfire_ID"] <- IPM_species_Medfire_ID[dom_species]
		} else {
			IPM_Medfire_Interface[i,"IPM_LCT_Medfire_ID"] <- 0
		}	
}
 	save(IPM_Medfire_Interface, file="./mdl_coupling/IPM_Medfire_Interface.rdata")
 }