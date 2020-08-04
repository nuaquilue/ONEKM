## 9. AFFORESTATION
      if(MEDFIRE){
        if(processes[afforest.id] & t %in% temp.afforest.schedule){
            aux  <- afforestation(land, coord, orography, clim, sdm)
            if (length(aux)!=0 & nrow(aux)>0){
              land$spp[land$cell.id %in% aux$cell.id] <- aux$spp
              land$age[land$cell.id %in% aux$cell.id] <- 0
              land$tsdist[land$cell.id %in% aux$cell.id] <- 0
              land$distype[land$cell.id %in% aux$cell.id] <- afforest
              clim$spp[clim$cell.id %in% aux$cell.id] <- aux$spp
              clim$sdm[clim$cell.id %in% aux$cell.id] <- 1
              clim$sqi[clim$cell.id %in% aux$cell.id] <- aux$sqi
              track.afforest <- rbind(track.afforest, data.frame(run=irun, year=t, table(aux$spp)))
            }
            temp.afforest.schedule <- temp.afforest.schedule[-1] 
        }
      }

      if(IPM) {
        if(COLONIZATION){
          # keep a dataframe with the basal area for potential colonization
          temp.suitable <- data.frame(ID = map$ID,
                                      saplings = integer(nrow(map)),
                                      plot_basal_area = integer(nrow(map)),
                                      #neigh = integer(nrow(map)),
                                      sp_age = integer(nrow(map)),
                                      neigh_ba = integer(nrow(map)),
                                      suitable = logical(nrow(map)))
          temp.suitable$plot_basal_area <- rowSums(ba) #sapply(X=trees,FUN=function(x) sum(unlist(x$ba)))
          # default 
          temp.suitable$suitable <- T
          colonized.plots.ID <- c()
          order <- sample(x = 1:NUM_SP,size = NUM_SP,replace = F)
          for(i in order){            
            print(paste(date()," - ",scn.name," - ",clim.scn," - year ",iyear," - colonization for sp ",i,"...",sep=""))
            temp.suitable$sp_age <- IPM.forest.age[,i]
            suitable <- GetSuitablePlots(map=map, temp.suitable=temp.suitable, ba = ba, distancias=distancias, sp=tesauro[i,1],
                                          max.dist = max.dist, verbose=T)   
            ##apply colonization
            if(!is.null(nrow(suitable))){
              if(nrow(suitable)>0){
              	if (MEDFIRE){
              		##Only suitable if it has been colonized by Medfire as well
              		if (processes[afforest.id]){
              			suitable <- left_join(suitable,map[,c("ID","Medfire.id")], by="ID") %>% 
              						left_join(land[,c("cell.id","spp")], by=c("Medfire.id"="cell.id"))
              			suitable <- suitable[(suitable$spp != 14),] #LCT that remain shrub after medfire colonization cannot be colonized
              		}
              	} ##if Medfire
                if(i %in% conifers) my.model <- colonization.glm[[1]]
                else if(i %in% quercus) my.model <- colonization.glm[[2]]
                else if(i %in% deciduous) my.model <- colonization.glm[[3]]
                suitable$colonized <- predict(my.model,newdata = suitable,type = "response")
                suitable$colonized <- ifelse(suitable$colonized > colonization.threshold,1,0)
                suitable<- suitable[suitable$colonized==1,]
                if (nrow(suitable)>0){
	                k <- match(suitable$ID,map$ID)
	                colonized.plots.ID <- c(colonized.plots.ID, map$ID[k])
	                #add new saplings to saplings list
	                for(j in 1:nrow(suitable)){
	                  saplings[k[j],tesauro[i,1]] <- new.saplings[tesauro[i,1]] 
	                }
            	}
              }
            }# if there are suitable plots
          }#for each species
          ##colonize plots colonized in Medfire but not in IPM
          colonized.plots.ID <- unique(colonized.plots.ID)
          if (MEDFIRE){
            if (processes[afforest.id]){
            	#colonized.medfire <- which(map$Medfire.id %in% aux$cell.id)
      			  #to.colonize.ID <- map$ID[colonized.medfire][!(map$ID[colonized.medfire] %in% colonized.plots.ID )] 
        			to.colonize.IPM.i <- which((map$Medfire.id %in% aux$cell.id) & !(map$ID %in% colonized.plots.ID )) 
        			if(length(to.colonize.IPM.i)>0){
        				#k <- match(to.colonize.ID,map$ID)
        				for(j in to.colonize.IPM.i){
        					spp.Medfire <- land$spp[land$cell.id==map$Medfire.id[j]]
        					spp.IPM <- Medfire.index.IPM.spp[spp.Medfire]
        					saplings[j, spp.IPM] <- new.saplings[spp.IPM]
        				}
        			}
      		  }
          } ## if Medfire
        }#if colonization module is active
      }#if IPM

      ## Track colonized plots LCT evolution of all the cells that have IPM plots
      if(MEDFIRE){
      	if(IPM.afforestation){
      		IPM.colonized.plots.indexes <- which(map$Medfire.id %in% land$cell.id[land$distype==afforest & !is.na(land$distype)])
        		for (col.plot in IPM.colonized.plots.indexes) {
        			if (sum(ba[col.plot,])>0){
        				IPM.spp <- which.max(ba[burnt.plot,])
        				Medfire.spp <- IPM.index.Medfire.spp[IPM.spp]
        				land$spp[land$cell.id %in%  map$Medfire.id[col.plot]] <- Medfire.spp
        				clim$spp[clim$cell.id %in%  map$Medfire.id[col.plot]] <- Medfire.spp
        			} else if(sum(saplings[col.plot,])>0){
        				IPM.spp <- which.max(saplings[burnt.plot,])
        				Medfire.spp <- IPM.index.Medfire.spp[IPM.spp]
        				land$spp[land$cell.id %in%  map$Medfire.id[col.plot]] <- Medfire.spp
        				clim$spp[clim$cell.id %in%  map$Medfire.id[col.plot]] <- Medfire.spp
        			} else{
        				land$spp[land$cell.id %in%  map$Medfire.id[col.plot]] <- 14
        				clim$spp[clim$cell.id %in%  map$Medfire.id[col.plot]] <- 14
        			}
        		}#for colonized plot
      	}
      }          	