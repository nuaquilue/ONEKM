burn.IPM <- function(){
	
	load(paste0("./outputs/",scn.name,"burnt_cells.rdata" ))
	load("./IPM/map_medfireID.rdata")
	load(adult.trees.file); load(ba.file); load(saplings.file)
	count = 0
	for (cell in burnt.cells[]){
		cell_index_IPM <- which(map.medfireID$cell_id_medfire == cell)
		if (length(cell_index_IPM)!=0){
			count= count +1
			for (s in 1:NUM_SP){
				adult.trees[[s]][cell_index_IPM,]<-0
				ba[cell_index_IPM,]<-0
				saplings[cell_index_IPM,]<-0
			} #for species 
		}#if burnt cell is in IPM map	
	}# for burnt cell
save(adult.trees, file=adult.trees.file)
save(ba, file=ba.file)
save(saplings, file=saplings.file)
}