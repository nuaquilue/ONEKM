##IPM-Medfire global variables
IPM <- T
MEDFIRE <- T

## Time lenght (in years) of a model simulation, from 2000 to 2100 
time.horizon <- 25

clim.scn <- "rcp45"
clim.mdl <- "SMHI-RCA4_MOHC-HadGEM2-ES"

nrun <- 1
testing <- F
# work.path <- "C:/Users/uriso/Desktop/ONEKM" #Laptop oriol
#work.path <- "C:/Users/Administrator/Oriol/ONEKM" #remote brotons lab
work.path <- "/home/xrdpuser/Desktop/ONEKM"
##Species interface vectors	                  
Medfire.index.IPM.spp<-c(5,6,8,9,7,10,2,12,15,11,11,3,1) #quercus humilis and faginea are classified as the same for IPM
IPM.index.Medfire.spp <-c(13,13,12,13,1,2,5,3,4,6,11,8,13,13,9,13) #classified as other (13): coniferes(1), decideous(2), juniperus(4), quercus pyrenaica(13), quercus robur(14) and 	Sclerophyllous (16)

out.path <- paste0("mdl_interface/output/", scn.name)
if(!file.exists(out.path)){
    dir.create(file.path(getwd(), out.path), showWarnings = T) 
}