rm(list=ls())
gc()
setwd("C:/Users/uriso/Desktop")

n.intervals.mesh <- 4000
study.area <- "BCN"
load(paste0("ONEKM/mdl_interface/IPM_",study.area, "_map.rdata"))
orig.adult.trees.file <- paste("ONEKM/IPM/initial_variables/trees_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
orig.ba.file <- paste("ONEKM/IPM/initial_variables/ba_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
orig.saplings.file <- paste("ONEKM/IPM/initial_variables/saplings_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
orig.plots.age.file <- paste("ONEKM/IPM/initial_variables/orig_plots_age", study.area, ".rdata",sep="")


load(orig.ba.file); load(orig.saplings.file)
load("ONEKM/Medfire/inputlyrs/rdata/land.rdata")
ini_ba <- ba
ini_sapl<- saplings
ini_land <- land

fire.regeneration <- c(T,T,T,F,T,F,T,F,F,F,T,T,T,T,T,T)
Medfire.index.IPM.spp<-c(5,6,8,9,7,10,2,12,15,11,11,3,1) #quercus humilis and faginea are classified as the same for IPM
IPM.index.Medfire.spp <-c(13,13,12,13,1,2,5,3,4,6,11,8,13,13,9,13)
scenario <- "coupled_scn_1_run_2"
##Test function to assess if for all burnt cells in Medfire within IPM study area
##The cells of IPM are also burnt
analyse.burnt.IPM.cells <- function(final.year){
	check <- 1
	for (iyear in  2000:final.year){
		hist_fires <- raster(paste0("ONEKM/Medfire/historic_fires/Fires_",iyear,".TIF"))
		burnt.cells <- which(!is.na(hist_fires[]))
		burnt.cells.bcn.IPM <- which(map$Medfire.id %in% burnt.cells)
		if (length(burnt.cells.bcn.IPM!=0)){
			ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
			load(ba.file)
			if(sum(ba[burnt.cells.bcn.IPM,])!=0){
				cat(paste0("Model does not burn IPM cell at year ", iyear, "\n"))
				check<-0
			}
		}
		cat(paste0("Analysing year: ", iyear, "\n"))
	}
	if (check) {
		cat("model can burn all IPM cells \n")
	}
}

##Test function to assess:
## 1. That after fire there are no saplings of species that cannot regenerate
## 2. That after fire there are saplings of species that can regenerate if they were present before (however, age )
analyse.IPM.saplings.post.fire<- function(final.year){
	check1 <-1 
	check2 <-1
	for (iyear in  2000:final.year){
	  cat(paste0("Analysing year: ", iyear, "\n"))
		hist_fires <- raster(paste0("ONEKM/Medfire/historic_fires/Fires_",iyear,".TIF"))
		burnt.cells <- which(!is.na(hist_fires[]))
		burnt.cells.bcn.IPM <- which(map$Medfire.id %in% burnt.cells)
		cat(paste0("Analysing ",length(burnt.cells.bcn.IPM), " burnt cells\n"))
		if (length(burnt.cells.bcn.IPM!=0)){
			##load sapling file
			saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear, "_run_1.rdata",sep="")
			load(saplings.file)
			##load previous ba file
			if (iyear==2000){
				ba.file <- paste("ONEKM/IPM/initial_variables/ba_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
			}
			else{
				ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear-1, "_run_1.rdata",sep="")
			}
			load(ba.file)
			##check previous species
			for (cell in burnt.cells.bcn.IPM){
			  #cat("running for loop \n")
				if (sum(ba[cell,])!=0){
					##check if there is presence of trees that can regenerate that were present before
					prev.trees.spp <- which(ba[cell,]>0)
					prev.trees.reg <- prev.trees.spp[prev.trees.spp %in% which(fire.regeneration==T)]
					if (!all(saplings[cell, prev.trees.reg]>0)){
						cat(paste0("Some trees that can regenerate don't appear in year ", iyear," and burnt cell: ", cell, "\n"))
						check1 <- 0
					}
					##check if there is presence of trees that cannot regenerate
					prev.trees.no.reg <- which(fire.regeneration==F)
					if(!all(saplings[cell, prev.trees.no.reg]==0)){
						cat(paste0("Some trees that cannot regenerate appear in year ", iyear," and burnt cell: ", cell, "\n"))
						check2 <- 0
					}
				} else if(!all(saplings[cell,]==0)){
					cat(paste0("Some trees appear in year ", iyear," and burnt cell: ", cell, "and no prev trees were present \n"))
						check2 <- 0
				}	
			}
		}
	}
	if(check1){
		cat("model predicts correctly for all plots that after fire there are saplings of species that can regenerate and were present before \n")
	}
	if(check2){
		cat("model works well since it does not put saplings that should not be there \n")
	}
}


##Medfire post fire regeneration just after fire, checks:
## 1: tsdist is updated to 0
## 2: tburnt is equal or bigger than 1
## 3: Medfire age for vegetations cells is 0
## 4: Spp after fire is the most abundant saplings (previous dominant tree species with regeneration traits)
analyse.Medfire.postfire <- function(final.year){
	check <- 1
	for (iyear in  2010:final.year){
	 	cat(paste0("Analysing year: ", iyear, "\n"))
		land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear, "_run_1.rdata",sep="")
		load(land.file)
		saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear, "_run_1.rdata",sep="")
		load(saplings.file)
		ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
		load(ba.file)

		hist_fires <- raster(paste0("ONEKM/Medfire/historic_fires/Fires_",iyear,".TIF"))
		burnt.cells <- which(!is.na(hist_fires[]))
		burnt.cells.bcn.IPM <- which(map$Medfire.id %in% burnt.cells)

		if (!all(land$tsdist[land$cell.id %in% burnt.cells]==1)){
			cat(paste0("Model does not update land tsdist ", iyear, "\n"))
		}
		if (!all(land$age[which(land$spp<=14 & (land$cell.id %in% burnt.cells))]==1)){
			cat(paste0("Model does not update land age ", iyear, "\n"))
		}
		if (!all(land$tburnt[land$spp<18 & (land$cell.id %in% burnt.cells)]>=1)){
			cat(paste0("Model does not update land tburnt ", iyear, "\n"))
		}

		if (length(burnt.cells.bcn.IPM!=0)){
			for (cell in burnt.cells.bcn.IPM){
			  if (sum(saplings[cell,])==0){
			    if (land$spp[land$cell.id %in% map$Medfire.id[cell]]!=14){
			      cat(paste0("Model does not correctly update land spp of burnt cell ", cell, " year ", iyear, "\n"))
			      check <-0
			      }
			  } else if (land$spp[land$cell.id %in% map$Medfire.id[cell]]!=IPM.index.Medfire.spp[which.max(saplings[cell,])]){
					cat(paste0("Model does not correctly update land spp of burnt cell ", cell, " year ", iyear, "\n"))
					check <-0
				}
			}
		}
		
	}
	if (check) {
		cat("model correctly updates land spp after fire \n")
	}
	
}

analyse.Medfire.dyn.spp <- function(final.year){
	check <- 1
	for (iyear in  2010:(final.year-10)){
		cat(paste0("Analysing year: ", iyear, "\n"))
		hist_fires <- raster(paste0("ONEKM/Medfire/historic_fires/Fires_",iyear,".TIF"))
		burnt.cells <- which(!is.na(hist_fires[]))
		burnt.cells.bcn.IPM <- which(map$Medfire.id %in% burnt.cells)
		if (length(burnt.cells.bcn.IPM!=0)){

			land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear, "_run_1.rdata",sep="")
			load(land.file)
			land0<-land
			saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear, "_run_1.rdata",sep="")
			load(saplings.file)
			saplings0 <- saplings
			ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
			load(ba.file)
			ba0 <- ba

			##Check that after 9 years saplings are the same as year 0 after fire and ba is still 0
			saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear + 9, "_run_1.rdata",sep="")
			load(saplings.file)
			ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear + 9, "_run_1.rdata",sep="")
			load(ba.file)
			land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear+10, "_run_1.rdata",sep="")
			load(land.file)
			burnt_once_IPM <- burnt.cells.bcn.IPM[map$Medfire.id[burnt.cells.bcn.IPM] %in% land$cell.id[land$tburnt==land0$tburnt]]
			if( !all(ba[burnt_once_IPM,]==0) | !all(saplings[burnt_once_IPM,]==saplings0[burnt_once_IPM,]) ){
				cat(paste0("Model changes saplings or ba after 9 years of land burnt in year", iyear, "\n"))
			}
			if (!all(land0$spp[land$cell.id %in% map$Medfire.id[burnt_once_IPM]] == land$spp[land$cell.id %in% map$Medfire.id[burnt_once_IPM]]))
				cat(paste0("Model changes land spp after 9 years of land burnt in year", iyear, "\n"))
			##Check that after 10 years if ba of trees is not 0 and saplings are the same as year 0 after fire
			saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear + 10, "_run_1.rdata",sep="")
			load(saplings.file)
			ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear + 10, "_run_1.rdata",sep="")
			load(ba.file)
			land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear+11, "_run_1.rdata",sep="")
			load(land.file)
			forest.cells <- burnt.cells.bcn.IPM[map$Medfire.id[burnt.cells.bcn.IPM] %in% land$cell.id[land0$spp<14 & land$tburnt>=1 & land$tburnt==land0$tburnt]]
			if( length(forest.cells)!= 0){
				cat(paste0("analysing ",length(forest.cells), " forest cells \n"))
				for (cell in forest.cells){
					if (sum(ba[cell,])==0 & sum(saplings0[cell,])!=0)
						cat(paste0("Model have BA of 0 in forest cells after 10 years of land burnt in year ", iyear, "\n"))
					## if check land spp equal ba
				  if (land$spp[land$cell.id %in% map$Medfire.id[cell]]!=IPM.index.Medfire.spp[which.max(ba[cell,])]){
				    cat(paste0("Model does not correctly update land spp of burnt cell ", cell, " year ", iyear, "\n"))
				    check <-0
				  }
				}
			}
			burnt_once_IPM <- burnt.cells.bcn.IPM[map$Medfire.id[burnt.cells.bcn.IPM] %in% land$cell.id[land$tburnt==land0$tburnt]]
			if (!all(saplings[burnt_once_IPM,]==saplings0[burnt_once_IPM,])){
				cat(paste0("Model changes saplings after 10 years of land burnt in year", iyear, "\n"))
			}
		}
	}
}
##IPM to Medfire
analyse.colonization <- function(final.year){
	for (iyear in  2010:final.year){
		col.file <- paste("remote_output/", scenario,"/colonized.plots.ID_", iyear, ".rdata",sep="")
	  load(col.file) 
	  land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear, "_run_1.rdata",sep="")
	  load(land.file)
	  cat(paste0("Year ", iyear, " - Analysing ", length(colonized.plots.ID), " colonized plots\n"))
	  col.IPM.index <- which(map$ID %in% colonized.plots.ID)
	  col.Medfire.cell.id <- land$cell.id[land$cell.id %in% map$Medfire.id[col.IPM.index]]
	  if (!all(land$tsdist[land$cell.id %in% col.Medfire.cell.id]==1)){
	    cat(paste0("Medfire does not update tsdist after colonization - year ", iyear))
	  }
	  if (!all(land$distype[land$cell.id %in% col.Medfire.cell.id]==10)){
	    cat(paste0("Medfire does not update disttype after colonization - year ", iyear))
	  }
	  
	  saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear - 1, "_run_1.rdata",sep="")
	  load(saplings.file)
	  ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear - 1, "_run_1.rdata",sep="")
	  load(ba.file)
	  new.col.cells <- which(apply(ba[col.IPM.index,],1,sum)==0 & apply(saplings[col.IPM.index,],1,sum)==0)
	  new.col.Medfire.cell.id <- land$cell.id[land$cell.id %in% map$Medfire.id[col.IPM.index[new.col.cells]]]
	  # if (!all(land$age[land$cell.id %in% new.col.Medfire.cell.id]==1)){
	  #   cat(paste0("Medfire does not update age after new colonized cells - year ", iyear))
	  # } ## only active if initial shurb IPM same as initial shrub Medfire
	  for (j in 0:floor((2089-iyear)/10)){
	    land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear+10*j, "_run_1.rdata",sep="")
	    load(land.file)
	    saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear+10*j, "_run_1.rdata",sep="")
	    load(saplings.file)
	    ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear+10*j, "_run_1.rdata",sep="")
	    load(ba.file)
	    curr_ba <- ba
	    ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear+10*j - 1, "_run_1.rdata",sep="")
	    load(ba.file)
	    prev_ba <- ba
	    ba_tot<-0
	    sapl_tot<-0
	    for (i in col.IPM.index){
	      if (sum(prev_ba[i,])!=0 & sum(curr_ba[i,])!=0){
	        ba_tot <- ba_tot +1 
	        if (land$spp[land$cell.id == map$Medfire.id[i]] != IPM.index.Medfire.spp[which.max(prev_ba[i,])]){
	          cat(paste0("Medfire does not track IPM spp cell with ba>0 ", i, " year ", j, " \n"))
	        }
	      }
	      else if (sum(saplings[i,])!=0){
	        sapl_tot <- sapl_tot +1
	        if(land$spp[land$cell.id == map$Medfire.id[i]] != IPM.index.Medfire.spp[which.max(saplings[i,])])
	          cat(paste0("Medfire does not track IPM spp cell with saplings>0 ", i, " year ", j, " \n"))
	      }
	      else if (land$spp[land$cell.id == map$Medfire.id[i]] != 14){
	        cat(paste0("Medfire does not track IPM spp empty cell ", i, " year ", j, " \n"))
	      }
	    }#for col cell
	    #print(paste0("ba tot: ", ba_tot, " sapl tot: ", sapl_tot))
	 }# for future year
  }# for current year
}

##Redo above functions for 90 years and colonization
##Make sure that colonized cells info pass to Medfire
##
ini_ba_col <- ini_ba[which(map$ID %in% colonized.plots.ID),]
length(which(apply(ini_ba_col,1,sum)==0))


##Medfire post fire regeneration just after fire, checks:
## 1: tsdist is updated to 0
## 2: tburnt is equal or bigger than 1
## 3: Medfire age for vegetations cells is 0
## 4: Spp after fire is the most abundant saplings (previous dominant tree species with regeneration traits)
analyse.Medfire.postfire.2 <- function(final.year){
	check <- 1
	for (iyear in  2010:final.year){
	  land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear-1, "_run_1.rdata",sep="")
	  load(land.file)
	  prev_land <- land
		land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear, "_run_1.rdata",sep="")
		load(land.file)
		saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear, "_run_1.rdata",sep="")
		load(saplings.file)
		ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
		load(ba.file)

		
		burnt.cells <- land$cell.id[(land$tburnt - prev_land$tburnt) == 1 & !is.na(land$tburnt - prev_land$tburnt)]
		burnt.cells.bcn.IPM <- which(map$Medfire.id %in% burnt.cells)

		cat(paste0("Analysing ", length(burnt.cells.bcn.IPM)," year: ", iyear, "\n"))
		
		if (!all(land$age[which(land$spp<=14 & (land$cell.id %in% burnt.cells) & land$distype==6)]==1)){
			cat(paste0("Model does not update land age ", iyear, "\n"))
		}
		
		for (j in 0:floor((2089-iyear)/10)){
		    land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear+10*j, "_run_1.rdata",sep="")
		    load(land.file)
		    saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear+10*j, "_run_1.rdata",sep="")
		    load(saplings.file)
		    ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear+10*j, "_run_1.rdata",sep="")
		    load(ba.file)
		    curr_ba <- ba
		    ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear+10*j - 1, "_run_1.rdata",sep="")
		    load(ba.file)
		    prev_ba <- ba
		    ba_tot<-0
		    sapl_tot<-0
		    for (i in burnt.cells.bcn.IPM){
		      if (prev_land$spp[land$cell.id %in% map$Medfire.id[i]])
		      if (sum(prev_ba[i,])!=0 & sum(curr_ba[i,])!=0){
		        ba_tot <- ba_tot +1 
		        if (land$spp[land$cell.id %in% map$Medfire.id[i]] != IPM.index.Medfire.spp[which.max(prev_ba[i,])]){
		          cat(paste0("Medfire does not track IPM spp cell with ba>0 ", i, " year ", j, " \n"))
		        }
		      }
		      else if (sum(saplings[i,])!=0){
		        sapl_tot <- sapl_tot +1
		        if(land$spp[land$cell.id %in% map$Medfire.id[i]] != IPM.index.Medfire.spp[which.max(saplings[i,])])
		          cat(paste0("Medfire does not track IPM spp cell with saplings>0 ", i, " year ", j, " \n"))
		      }
		      else if (land$spp[land$cell.id %in% map$Medfire.id[i]] != 14){
		        cat(paste0("Medfire does not track IPM spp empty cell ", i, " year ", j, " \n"))
		      }
		    }#for burnt cell
		    #print(paste0("ba tot: ", ba_tot, " sapl tot: ", sapl_tot))
		}# for future year
	}#for current year
}