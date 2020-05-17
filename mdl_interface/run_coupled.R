IPM.Medfire.mdl <- function(scn.name){
	source("./mdl_coupling/init_IPM.r")
	source("./mdl_coupling/init_Medfire.r")
	source("./mdl_coupling/Medfire_step.r")
  	source("./mdl_coupling/IPM_step.r")
	## Start the simulations   
	irun=1   # for testing
	for(irun in 1:nrun){
    
	    ## Copy the schedulings in auxiliar vectors (only for those processes included in the current version)
	    temp.clim.schedule <- clim.schedule
	    temp.lchg.schedule <- lchg.schedule
	    temp.fmgmt.schedule <- fmgmt.schedule
	    temp.fire.schedule <- fire.schedule
	    temp.pb.schedule <- pb.schedule
	    temp.drought.schedule <- drought.schedule
	    temp.post.fire.schedule <- post.fire.schedule
	    temp.cohort.schedule <- cohort.schedule
	    temp.afforest.schedule <- afforest.schedule
	    temp.growth.schedule <- growth.schedule
	    
	    
	    ## Load initial spatial dynamic state variables in a data.frame format
	    load("inputlyrs/rdata/land.rdata")
	    save(land, file="./mdl_coupling/output/land.rdata")
	    iyear= 2000
	    for(t in time.seq){
	    	IPM.step(iyear)
	    	Medfire.step(t)
	    	burn_IPM()
	    	iyear = iyear+1

	    }
		## Print maps at the end of the simulation period per each run
	    # if(write.sp.outputs){
	    #   BURNT[!is.na(MASK[])] <- land$tburnt
	    #   writeRaster(BURNT, paste0(out.path, "/lyr/TimesBurnt_r", irun, ".tif"), format="GTiff", overwrite=T)
	    # }      

	    ######include if all results IPM##############
	} #run

cat("... writing outputs", "\n")
write.table(track.fmgmt[-1,], paste0(out.path, "/Management.txt"), quote=F, row.names=F, sep="\t")
track.fire$rem <- pmax(0,track.fire$atarget-track.fire$aburnt.highintens-track.fire$aburnt.lowintens)
write.table(track.fire[-1,], paste0(out.path, "/Fires.txt"), quote=F, row.names=F, sep="\t")
write.table(track.pb[-1,], paste0(out.path, "/PrescribedBurns.txt"), quote=F, row.names=F, sep="\t")
write.table(track.drougth[-1,], paste0(out.path, "/Drought.txt"), quote=F, row.names=F, sep="\t")
names(track.post.fire)[4:5] <- c("spp.in", "ha")
write.table(track.post.fire[-1,], paste0(out.path, "/PostFire.txt"), quote=F, row.names=F, sep="\t")
names(track.cohort)[4:5] <- c("spp.in", "ha")
write.table(track.cohort[-1,], paste0(out.path, "/Cohort.txt"), quote=F, row.names=F, sep="\t")
names(track.afforest)[3:4] <- c("spp", "ha")
write.table(track.afforest[-1,], paste0(out.path, "/Afforestation.txt"), quote=F, row.names=F, sep="\t")
write.table(track.land[-1,], paste0(out.path, "/Land.txt"), quote=F, row.names=F, sep="\t")
}