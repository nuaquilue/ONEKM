##########
##  READ CLIMATE IPM
## Output: dataframe Medfire.id | ID | X_UTM | Y_UTM | temp | precip | anom_temp | anom_pl
## Anom temp: T(study period)-T(reference period)
## anom_pl: (P(study period)-P(reference period))/P(reference period)
##########

read.IPM.climate <- function(work.path){
  
	## Mask of the new study area
	load(paste0(work.path,"/inputlyrs/rdata/mask.rdata"))
	#MASK <- MASK_BCN
	load(paste0(work.path,"/mdl_interface/cell.id.interface.rdata"))
	extCat <- extent(c(250000, 540000, 4480000, 4760000))

	for(clim.scn in c("rcp45", "rcp85")){
		for(clim.mdl in c("KNMI-RACMO22E_ICHEC-EC-EARTH",
	                      "KNMI-RACMO22E_MOHC-HadGEM2-ES",
	                      "SMHI-RCA4_CNRM-CERFACS-CNRM-CM5",
	                      "SMHI-RCA4_MPI-M-MPI-ESM-LR",
	                      "SMHI-RCA4_MOHC-HadGEM2-ES")){

			clima <- data.frame(Medfire.id=1:ncell(MASK), Mask=MASK[])

			decade <- "Hist19712000"
			var <- "TNMM"   ##TXMM should be mean temperature
			TEMP.hist.min <- raster(paste0(work.path, "/inputlyrs/asc/ClimDownscaled/", var,"_",  clim.mdl, "_", decade, "_1000m.asc"))
			TEMP.hist.min <- extend(TEMP.hist.min, extCat)
			var <- "TXMM"   ##TXMM should be mean temperature
			TEMP.hist.max <- raster(paste0(work.path, "/inputlyrs/asc/ClimDownscaled/", var,"_",  clim.mdl, "_", decade, "_1000m.asc"))
			TEMP.hist.max <- extend(TEMP.hist.max, extCat)
			TEMP.hist <- (TEMP.hist.max[] + TEMP.hist.min[])/2

			var <- "PRCPTOT" ##MEAN ANNUAL TEMPERATURE
			PRECIP.hist <- raster(paste0(work.path, "/inputlyrs/asc/ClimDownscaled/", var,"_",  clim.mdl, "_", decade, "_1000m.asc"))
			PRECIP.hist <- extend(PRECIP.hist, extCat)


			##first constructs decade 00
			var <- "TNMM"
			TEMP.min <- raster(paste0(work.path, "/inputlyrs/asc/ClimDownscaled/", var,"_",  clim.scn, "_", clim.mdl,"_proj00_1000m.asc"))
			TEMP.min <- extend(TEMP.min, extCat)
			var <- "TXMM"
			TEMP.max <- raster(paste0(work.path, "/inputlyrs/asc/ClimDownscaled/", var,"_",  clim.scn, "_", clim.mdl,"_proj00_1000m.asc"))
			TEMP.max <- extend(TEMP.max, extCat)
			TEMP <- (TEMP.max[] + TEMP.min[])/2
			var <- "PRCPTOT"
			PRECIP <- raster(paste0(work.path, "/inputlyrs/asc/ClimDownscaled/", var,"_",  clim.scn, "_", clim.mdl,"_proj00_1000m.asc"))
			PRECIP <- extend(PRECIP, extCat)
			clima$temp <- TEMP[]
			clima$precip <- PRECIP[]
			clima$anom_temp <- TEMP[] - TEMP.hist[]
			clima$anom_pl <- (PRECIP[] - PRECIP.hist[])/PRECIP.hist[]

			clima <- filter(clima, Mask==1) %>% select(-Mask) %>% left_join(cell.id.interface[,c("Medfire.id","IFN.id")], by="Medfire.id")
			#clima <- clima[!is.na(clima$IFN.id),]
			names(clima)[names(clima) == 'IFN.id'] <- 'ID'

			decade <- "00"
			clima_IPM <- clima
			save(clima_IPM, file=paste0(work.path, "/mdl_interface/inputlyrs/rdata/climate_", clim.scn, "_", clim.mdl, "_", decade, ".rdata"))

			for(decade in seq(10,90,10)){ 
				clima <- data.frame(Medfire.id=1:ncell(MASK), Mask=MASK[])

				print(paste("Building: scenario", clim.scn, "model", clim.mdl, "- decade", decade))


				var <- "TNMM"
				TEMP.min <- raster(paste0(work.path, "/inputlyrs/asc/ClimDownscaled/", var,"_",  clim.scn, "_", clim.mdl, "_proj", decade, "_1000m.asc"))
				TEMP.min <- extend(TEMP.min, extCat)
				var <- "TXMM"
				TEMP.max <- raster(paste0(work.path, "/inputlyrs/asc/ClimDownscaled/", var,"_",  clim.scn, "_", clim.mdl, "_proj", decade, "_1000m.asc"))
				TEMP.max <- extend(TEMP.max, extCat)
				TEMP <- (TEMP.max[] + TEMP.min[])/2
				var <- "PRCPTOT"
				PRECIP <- raster(paste0(work.path, "/inputlyrs/asc/ClimDownscaled/", var,"_",  clim.scn, "_", clim.mdl, "_proj", decade, "_1000m.asc"))
				PRECIP <- extend(PRECIP, extCat)
				clima$temp <- TEMP[]
				clima$precip <- PRECIP[]
				clima$anom_temp <- TEMP[] - TEMP.hist[]
				clima$anom_pl <- (PRECIP[] - PRECIP.hist[])/PRECIP.hist[]

				clima <- filter(clima, Mask==1) %>% select(-Mask) %>% left_join(cell.id.interface[,c("Medfire.id","IFN.id")], by="Medfire.id")
				#clima <- clima[!is.na(clima$IFN.id),]
				names(clima)[names(clima) == 'IFN.id'] <- 'ID'
				clima_IPM <- clima
				save(clima_IPM, file=paste0(work.path, "/mdl_interface/inputlyrs/rdata/climate_", clim.scn, "_", clim.mdl, "_", decade, ".rdata"))

      }
    }
  }
  
}