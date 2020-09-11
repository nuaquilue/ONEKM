###############################################################################
#########single plot analysis#################################################
##############################################################################
setwd("C:/Users/uriso/Desktop/ONEKM")
source("./IPM/IPM_functions_v20_old.R")
library(Rcpp)
sourceCpp("./IPM/IPM_functions_v22_1_Year.cpp")

#load global and

trees.count <- function(trees, h, nx) {
  count <- rep(0, nrow(adult.trees[[1]]))
  for (i in 1:16){
    for(plot in 1:nrow(adult.trees[[1]])){
      count[plot] <- count[plot] + quadTrapezCpp_1(adult.trees[[i]][plot,],h[i],nx)
    }
  }
  return(count)
}
#dbh_dist <- matrix(nrow=92,ncol=nx)
#initial data
study.area <- "BCN"
n.intervals.mesh <- 2000
min.DBH <- 7.5  # in cm.

x2000 <- x.per.species(min.dbh=min.DBH,n.intervals=2000)
h2000 <- x2000[2,]-x2000[1,]
nx2000 <- 2000+1

x3000 <- x.per.species(min.dbh=min.DBH,n.intervals=3000)
h3000 <- x3000[2,]-x3000[1,]
nx3000 <- 3000+1

x4000 <- x.per.species(min.dbh=min.DBH,n.intervals=4000)
h4000 <- x4000[2,]-x4000[1,]
nx4000 <- 4000+1

orig.adult.trees.file <- paste("../ONEKM/IPM/initial_variables/trees_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
orig.ba.file <- paste("../ONEKM/IPM/initial_variables/ba_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
orig.saplings.file <- paste("../ONEKM/IPM/initial_variables/saplings_", study.area,"_", n.intervals.mesh, ".rdata",sep="")

load(orig.adult.trees.file);load(orig.ba.file); load(orig.saplings.file) ##make sure that it updates values

ini_cells <- c()
NUM_PLOTS <- 100
i <- 1
while(length(ini_cells)<NUM_PLOTS){
  if (sum(ba[i,])>0){
    ini_cells <- c(ini_cells, i)
  }
  i <- i +1
}
#map <- map[ini_cells,]
ba <- ba[ini_cells,]
saplings <- saplings[ini_cells,]
#IPM.forest.age <- IPM.forest.age[ini_cells,]
adult.trees<- lapply(adult.trees, function(x) {x[ini_cells,]})

ini_ba <- ba
ini_sapl <- saplings
ini_trees <- adult.trees


iyears <- seq(2000,2089,1)
ba_df <- data.frame(matrix(ncol = length(iyears)+1, nrow = 8))
#sapl_df <- data.frame(matrix(ncol = length(iyears)+1, nrow = 2))
rownames(ba_df) <- c("ba.1y.2000m", "ba.1y.3000m", "ba.1y.4000m","ba.10y", "num.trees.1y.2000m","num.trees.1y.3000m","num.trees.1y.4000m", "num.trees.10y")
colnames(ba_df) <- c("IFN3",iyears) 
#rownames(sapl_df) <- c("sapl.1y", "sapl.10y")
#colnames(sapl_df) <- c("IFN3",iyears) 
ba_df$IFN3[c(1,2,3,4)]<- apply(ini_ba,1,sum)[1]
ba_df$IFN3[c(5,6,7,8)]<- trees.count(ini_trees, h2000, nx2000)[1]
#sapl_df$IFN3[c(1,2)]<- apply(saplings,1,sum)[1]
#IPM data
setwd("C:/Users/uriso/Desktop/IPM_results")
scn.names <-c("only_growth_2000", "only_growth_3000", "only_growth_4000") #c("growth_death_ingrowth_2000", "growth_death_ingrowth_3000", "growth_death_ingrowth_4000")#c("growth_death_2000", "growth_death_3000", "growth_death_4000")# 
SIM_NUMBER <- c("only_growth_10y_2000")#c("growth_death_ingrowth_10y_2000")#c("growth_death_10y_2000")#
irun <- 1
ipm.10y.path <- "C:/Users/uriso/Desktop/IPM_10y" #oriol laptop
#for(sim in 1:length(versions)){ #ba_2000_1plot_only_growth

for (iyear in iyears){
  index <- 1
  for(scn.name in scn.names){
    adult.trees.file <- paste0("./", scn.name, "/trees_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
    ba.file <- paste0("./", scn.name, "/ba_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
    #saplings.file <- paste0("./saplings_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
    load(adult.trees.file);load(ba.file)# ;load(saplings.file) ##make sure that it updates values
    ba_df[index, as.character(iyear)]<- apply(ba,1,sum)[1]
    if (index==1) numb_trees<- trees.count(adult.trees, h2000, nx2000)
    else if (index==2 ) numb_trees<- trees.count(adult.trees, h3000, nx3000)
    else numb_trees<- trees.count(adult.trees, h4000, nx4000)
    ba_df[index+4, as.character(iyear)]<-numb_trees[1]
    #sapl_df[1, as.character(iyear)]<- apply(saplings,1,sum)[1]
    index <- index + 1
  }
  if ((iyear+1)%%10 == 0) {
    adult.trees.file <- paste("./", SIM_NUMBER,"/trees_",(iyear-9),"_",SIM_NUMBER, ".rdata",sep="")
    ba.file <- paste("./", SIM_NUMBER,"/ba_",(iyear-9),"_",SIM_NUMBER, ".rdata",sep="")
    #saplings.file <- paste(ipm.10y.path,"/resultados/saplings_",(iyear-9),"_",SIM_NUMBER, ".rdata",sep="")
    load(adult.trees.file);load(ba.file)#; load(saplings.file) ##make sure that it updates values
    ba_df[4, as.character(iyear) ]<- apply(ba,1,sum)[1]
    #sapl_df[8, as.character(iyear)]<- apply(saplings,1,sum)[1]
    ba_df[8, as.character(iyear) ]<- trees.count(adult.trees, h2000, nx2000)[1]
  }
}
title<-"Growth and survival: number of trees evolution"
##plot number of trees
plot(2000:2090, ba_df[7,], ylim=c(0, max(ba_df[5:8,]+5, na.rm=T)), type="l", lwd=2, col="red",lty=2, xlab="years", ylab="total number of trees")
lines(2000:2090,ba_df[5,], col="yellow", lwd=2, lty=2)
lines(2000:2090,ba_df[6,], col="green", lwd=2, lty=2)
lines(seq(2000,2090,10),ba_df[8, seq(1,100,10)], col="black", lwd=2)
legend("bottomleft", legend=c("1y_m_4000", "1y_m_3000", "1y_m_2000", "10y_m_2000"),
       col=c("red", "green", "yellow", "black"), lty=c(2,2,2,1), cex=0.8)
##only 2000
plot(2000:2090, ba_df[5,], ylim=c(0, max(ba_df[5:8,]+5, na.rm=T)), type="l", lwd=3, col="yellow",lty=2, xlab="years", ylab="total number of trees")
lines(seq(2000,2090,10),ba_df[8, seq(1,100,10)], col="black", lwd=3)
legend("bottomleft", legend=c("1y_m_2000", "10y_m_2000"),
       col=c("yellow", "black"), lty=c(2,1), cex=0.8)
#
##plot ba
title<-"Growth and survival: total basal area evolution"
plot(2000:2090, ba_df[3,], ylim=c(0, max(ba_df[1:4,]+5, na.rm=T)), type="l", lwd=2, col="red",lty=2, xlab="years", ylab="total basal area")
lines(2000:2090,ba_df[1,], col="yellow", lwd=2, lty=2)
lines(2000:2090,ba_df[2,], col="green", lwd=2, lty=2)
lines(seq(2000,2090,10),ba_df[4, seq(1,100,10)], col="black", lwd=2)
legend("bottomright", legend=c("1y_m_4000", "1y_m_3000", "1y_m_2000", "10y_m_2000"),
       col=c("red", "green", "yellow", "black"), lty=c(2,2,2,1), cex=0.8)
##only 2000
plot(2000:2090, ba_df[1,], ylim=c(0, max(ba_df[1:4,]+5, na.rm=T)), type="l", lwd=3, col="yellow",lty=2, xlab="years", ylab="total basal area")
lines(seq(2000,2090,10),ba_df[4, seq(1,100,10)], col="black", lwd=3)
legend("bottomright", legend=c( "1y_m_2000", "10y_m_2000"),
       col=c("yellow", "black"), lty=c(2,1), cex=0.8)
save(ba_df,file="./mesh.comparison.plot1.rdata")


#################
##mesh transform
#################

##2000 to 4000
iyear<-2089
ba.file <- paste("./", SIM_NUMBER,"/ba_",(iyear-9),"_",SIM_NUMBER, ".rdata",sep="")
adult.trees.file <- paste("./", SIM_NUMBER,"/trees_",(iyear-9),"_",SIM_NUMBER, ".rdata",sep="")
load(adult.trees.file);load(ba.file)#; load(saplings.file) ##make sure that it updates values
title<-"DBH distribution of pinus halepensis after 90 years of simulation"
plot(seq(1,4001,2), main = title, adult.trees[[5]][1,], ylim=c(0, max(adult.trees[[5]][1,])+1),
     lwd=2, type="l", col="yellow", xlim=c(0,2000))

colors<-c("yellow", "green", "red")
steps <-c(2,4/3,1)
i<-1
for(scn.name in scn.names){
  adult.trees.file <- paste0("./", scn.name, "/trees_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
  load(adult.trees.file)
  lines(seq(1,4001,steps[i]),adult.trees[[5]][1,], col=colors[i], lty=2)
  i<-i+1
}

legend("topright", legend=c("1y_m_4000", "1y_m_3000", "1y_m_2000", "10y_m_2000"),
       col=c("red", "green", "yellow", "yellow"), lty=c(2,2,2,1), cex=0.8)
##unncessesary code
# expanded_2000 <- list()
# for(i in 1:16){
# expanded_2000[[i]] <- matrix(0,1, nx4000)
# }
# 
# for (sp in 1:16){
# 	if ( sum(adult.trees[[sp]][1,])> 0 ){
# 		for (i in 1:2001){
# 			expanded_2000[2*i-1]<- adult.trees[[sp]][1,i]
# 			if(i!=2001){
# 			  expanded_2000[2*i]<- (adult.trees[[sp]][1,i] + adult.trees[[sp]][1,i+1])/2
# 			}
# 		}
# 	}
# }
# 
# ##3000 to 4000
# for(i in 1:16){
# expanded_3000[[i]] <- matrix(0,1, nx4000)
# }
# 
# for (sp in 1:16){
# 	if ( sum(adult.trees[[sp]][1,])> 0 ){
# 		diff<-0
# 		for (i in seq(1,3001,3) ){
# 			expanded_3000[i+diff]<- adult.trees[[sp]][1,i]
# 			if (i!=3001){
# 				slope <- (adult.trees[[sp]][1,i+1] - adult.trees[[sp]][1,i])
# 				expanded_3000[i+diff+1]<- slope*0.75+ adult.trees[[sp]][1,i]
# 				expanded_3000[i+diff+2]<- (adult.trees[[sp]][1,i+1] + adult.trees[[sp]][1,i+2])/2
# 				slope <- (adult.trees[[sp]][1,i+3] - adult.trees[[sp]][1,i+2])
# 				expanded_3000[i+diff+3]<- slope*0.25 + adult.trees[[sp]][1,i+2]
# 			}
# 			diff<-diff+1
# 		}
# 	}
# }