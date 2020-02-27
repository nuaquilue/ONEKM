######################################################################################
##
######################################################################################

drought <- function(land, clim, t){
  
  ## Tracking
  cat("Drought", "\n")
  
  ## Count how many ha to kill, totally and this time step
  to.kill <- filter(cbind(land, select(clim, sdm)), spp<=13, sdm==0, distype==0) 
  nkill <- round(table(to.kill$spp) /(10 - (t-(t %/% 10)*10) + 1))
  
  ## Kill randomly as many cells per spp
  killed.cells <- integer()
  for(i in names(nkill)){
    if(nkill[i]>0)
      kill.cells <- c(kill.cells, sample(to.kill$cell.id[to.kill$spp==i], nkill[i], replace=F))
  }
  
  return(killed.cells)
}
