######################################################################################
###  Land cover changes
###
######################################################################################

land.cover.change <- function(land, coord, orography, lc.trans=1, chg.cells){
                        
  cat(paste0("Land-cover transition: ", ifelse(lc.trans==1, "Urbanization.", ifelse(lc.trans==2, "Agriculture conversion.", 
      ifelse(lc.trans==3, "Rural abandonment.", "Undefined.")))))
  
  
  return(chg.cells)
}