
prob.igni <- function(land, orography, clim, interface){
  
  ## Tracking
  # print("Prob. ignition")
  
  ## Read ignition probability model
  pigni.mdl <- unlist(read.table("inputfiles/ProbIgniMdl.txt", header=T))
  
  z = pigni.mdl[1] + pigni.mdl[2]*orography$elev + pigni.mdl[3]*orography$slope/10 + pigni.mdl[4]*clim$precip + 
      pigni.mdl[5]*(interface$x == 3) + pigni.mdl[6]*(interface$x == 6)+ pigni.mdl[7]*(interface$x == 7) + 
      pigni.mdl[8]*orography$road/100  ## comprobar que slope/10 i road/100 estan ben dividits, unitats!!
  
  return((1/(1+exp(-1*z)))*100)
}