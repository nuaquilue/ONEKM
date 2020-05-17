######################################################################################################
## Reads burnt.cells by Medfire and using the cell.id interface finds the plots to be burnt in IPM
## For each burnt plot:
## 		- Calculates future.saplings to be recruited after 10 years.
## 		- Sets adult.trees, ba and saplings equal 0.
##		- Sets IPM.forest.age to 0.
######################################################################################################
## TO DO:
## Think what to do with low intensity burnt cells
######################################################################################################

burn.IPM <- function(){
	
	load(paste0("./outputs/",scn.name,"burnt_cells.rdata" ))
	load("./mdl_interface/cell.id.interface.rdata")

	##IPM plots burnt by medfire
	burnt.cells.IPM.index <- cell.id.interface$IPM.index[cell.id.interface$Medfire.id %in% burnt.cells]
	
	if(length(burnt.cells.IPM.index)!=0){
	
		load(adult.trees.file); load(ba.file); load(saplings.file); load(future.saplings.file); load(IPM.forest.age.file)
		##Calculate future saplings
		for (i in burnt.cells.IPM.index){
			##for debugging
			if (sum(ba[i,])==0){ ##When young forest or bush are burnt
				future.saplings[i,] <- 0	
			}
			else{
				if (IPM.forest.age[i]>9){
					for(j in 1:NUM_SP){
						if (fire.regeneration[j]){
							## future saplings: fixed number of new.saplings times the abundance proportion of the species in the plot (in ba) 
							future.saplings[i,j] <- new.saplings[j]*(ba[i,j]/sum(ba[i,]))
						} #if species can regenerate
						else {
							future.saplings[i,j] <- 0
						} #species cannot regenerate
					} #for tree species
				} #if plot age is older than 9
				else{
					future.saplings[i,] <- 0
				}
			}#if ba plot greater than 0
		} #for burnt IPM plot
		
		#Burns IPM plots by assigning adult.trees, ba and saplings to 0
		for (s in 1:NUM_SP){
				adult.trees[[s]][burnt.cells.IPM.index,]<-0
			} #for species 
		ba[burnt.cells.IPM.index,]<-0
		saplings[burnt.cells.IPM.index,]<-0
		IPM.forest.age[burnt.cells.IPM.index]<-0

		save(adult.trees, file=adult.trees.file)
		save(ba, file=ba.file)
		save(saplings, file=saplings.file)
		save(future.saplings, file=future.saplings.file)
		save(IPM.forest.age, file=IPM.forest.age.file)
	} #if burnt plots

}