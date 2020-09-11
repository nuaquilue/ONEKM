###############################################################################
#########single plot analysis#################################################
##############################################################################
rm(list=ls())
gc()
setwd("C:/Users/uriso/Desktop/ONEKM")
source("./IPM/IPM_functions_v20_old.R")
library(Rcpp)
library(ggplot2)
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

x5000 <- x.per.species(min.dbh=min.DBH,n.intervals=5000)
h5000 <- x5000[2,]-x5000[1,]
nx5000 <- 5000+1

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
ba_df <- data.frame(matrix(ncol = length(iyears)+1, nrow = 10))
#sapl_df <- data.frame(matrix(ncol = length(iyears)+1, nrow = 2))
rownames(ba_df) <- c("ba.1y.2000m", "ba.1y.3000m", "ba.1y.4000m", "ba.1y.5000m","ba.10y",
 "num.trees.1y.2000m","num.trees.1y.3000m","num.trees.1y.4000m", "num.trees.1y.5000m", "num.trees.10y")
colnames(ba_df) <- c("IFN3",iyears) 
#rownames(sapl_df) <- c("sapl.1y", "sapl.10y")
#colnames(sapl_df) <- c("IFN3",iyears) 
ba_df$IFN3[c(1,2,3,4,5)]<- mean(apply(ini_ba,1,sum))
ba_df$IFN3[c(6,7,8,9,10)]<- mean(trees.count(ini_trees, h2000, nx2000))
#sapl_df$IFN3[c(1,2)]<- apply(saplings,1,sum)[1]
#IPM data
setwd("C:/Users/uriso/Desktop/output")
scn.names <- c("g.d.i_100_plots_2000", "g.d.i_100_plots_3000", "100_plots_v2", "g.d.i_100_plots_5000")#c("only_growth_2000", "only_growth_3000", "only_growth_4000")#c("growth_death_2000", "growth_death_3000", "growth_death_4000")#
SIM_NUMBER <- c("IPM_10y_BCN_100_v3")#c("only_growth_10y_2000")#c("growth_death_10y_2000")
irun <- 1
ipm.10y.path <- "C:/Users/uriso/Desktop/IPM_10y/resultados/10y_v3" #oriol laptop
#for(sim in 1:length(versions)){ #ba_2000_1plot_only_growth

for (iyear in iyears){
  cat(iyear)
  index <- 1
  for(scn.name in scn.names){
    adult.trees.file <- paste0("./", scn.name, "/trees_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
    ba.file <- paste0("./", scn.name, "/ba_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
    #saplings.file <- paste0("./saplings_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
    load(adult.trees.file);load(ba.file)# ;load(saplings.file) ##make sure that it updates values
    ba_df[index, as.character(iyear)]<- mean(apply(ba,1,sum))
    if (index==1) numb_trees<- mean(trees.count(adult.trees, h2000, nx2000))
    else if (index==2 ) numb_trees<- mean(trees.count(adult.trees, h3000, nx3000))
    else if (index==3) numb_trees<- mean(trees.count(adult.trees, h4000, nx4000))
    else numb_trees<- mean(trees.count(adult.trees, h5000, nx5000))
    ba_df[index+5, as.character(iyear)]<-numb_trees[1]
    #sapl_df[1, as.character(iyear)]<- apply(saplings,1,sum)[1]
    index <- index + 1
  }
  if ((iyear+1)%%10 == 0) {
    adult.trees.file <- paste(ipm.10y.path,"/trees_",(iyear-9),"_",SIM_NUMBER, ".rdata",sep="")
    ba.file <- paste(ipm.10y.path,"/ba_",(iyear-9),"_",SIM_NUMBER, ".rdata",sep="")
    #saplings.file <- paste(ipm.10y.path,"/saplings_",(iyear-9),"_",SIM_NUMBER, ".rdata",sep="")
    load(adult.trees.file);load(ba.file)#; load(saplings.file) ##make sure that it updates values
    ba_df[5, as.character(iyear) ]<- mean(apply(ba,1,sum))
    #sapl_df[8, as.character(iyear)]<- apply(saplings,1,sum)[1]
    ba_df[10, as.character(iyear) ]<- mean(trees.count(adult.trees, h4000, nx4000))
  }
}
title<-"IPM: mean number of trees 100 plots" #main=title
##plot number of trees
colors<-c("green", "yellow", "orange", "red")
plot(seq(2000,2090,10), ba_df[10,seq(1,100,10)], ylim=c(0, max(ba_df[6:10,]+5, na.rm=T)), type="l", lwd=2, col="black", xlab="years", ylab="mean number of trees per plot")
for(i in 1:5){
  lines(2000:2090,ba_df[i+5,], col=colors[i], lwd=2, lty=2)
}
legend("bottomleft", legend=c("10y","1y_m_2000","1y_m_3000", "1y_m_4000", "1y_m_5000"),
       col=c("black", "green", "yellow", "orange","red"), lty=c(1,2,2,2,2,2), cex=0.8)

##plot ba
title<-"IPM: mean basal area 100 plots"
plot(seq(2000,2090,10), ba_df[5,seq(1,100,10)], ylim=c(0, max(ba_df[1:5,]+5, na.rm=T)), type="l", lwd=2, col="black", xlab="years", ylab="mean total BA per plot")
for(i in 1:5){
  lines(2000:2090,ba_df[i,], col=colors[i], lwd=2, lty=2)
}
legend("bottomright", legend=c("10y","1y_m_2000","1y_m_3000", "1y_m_4000", "1y_m_5000"),
       col=c("black", "green", "yellow", "orange","red"), lty=c(1,2,2,2,2,2), cex=0.8)



save(ba_df,file="./mesh.comparison.100plots.rdata")

##plot hist
rmse <- function(m, o){
  sqrt(mean((m - o)^2))
}

iyear<-2089
adult.trees.file <- paste(ipm.10y.path,"/trees_",(iyear-9),"_",SIM_NUMBER, ".rdata",sep="")
load(adult.trees.file)
final_nt_10y <- trees.count(adult.trees, h4000, nx4000)#trees.count(adult.trees)
index<-1
len_mesh <- c("m2000", "m3000", "m4000", "m5000")
rmse_vector<-numeric(4)
diff_df <- data.frame(diff=double(), mesh=integer())
for(scn.name in scn.names){
  adult.trees.file <- paste0("./", scn.name, "/trees_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
  load(adult.trees.file)
  if (index==1) numb_trees<- trees.count(adult.trees, h2000, nx2000)
  else if (index==2 ) numb_trees<- trees.count(adult.trees, h3000, nx3000)
  else if (index==3) numb_trees<- trees.count(adult.trees, h4000, nx4000)
  else numb_trees<- trees.count(adult.trees, h5000, nx5000)
  final_nt_1y <- numb_trees
  #[c(9,10,25,30,77,79,81,90,93,94,96)]
  #nt_diff <- data.frame(diff=(final_nt_1y/final_nt_10y*100-100))#(final_nt_10y - final_nt_1y)
  nt_diff <- data.frame(diff=(final_nt_1y-final_nt_10y))
  rmse_vector[index]<- rmse(final_nt_10y, final_nt_1y)
  nt_diff$mesh <- len_mesh[index]
  diff_df <- rbind(diff_df, nt_diff)
  #hist(nt_diff, breaks=100, main=paste0("mun trees diff at year ", iyear+1, scn.name))
  index<-index+1
}

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#
legend("topleft", legend = len_mesh, pch=15, pt.cex=3, cex=1.5, bty='n',
       col = colors)
##mean, median, variance, +-5, +-10, +-20
summary_df<- data.frame(mean=rep(0,4), median= rep(0,4), variance= rep(0,4),
                        a5 = rep(0,4), a10 = rep(0,4), a20= rep(0,4))
index <- 1
for (mesh in len_mesh){
  diff <- diff_df[which(diff_df$mesh==mesh),1]
  summary_df$mean[index]<- round(mean(diff),2 ) 
  summary_df$median[index]<- round(median(diff),2)
  summary_df$variance[index]<- round(var(diff),2)
  summary_df$a5[index]<- length(which(diff<5 & diff>-5))
  summary_df$a10[index]<- length(which(diff<10 & diff>-10))
  summary_df$a20[index]<- length(which(diff<20 & diff>-20))
  index <- index+1
}
diff_df$col<-1
i<-1
for (mesh in len_mesh){
  diff_df$col[which(diff_df$mesh==mesh)]<-colors[i]
  i<-i+1
}
write.table(summary_df, file="./summary_num_trees_no_percent.csv", dec=".", sep=";")
length(which(diff_df$diff<(-400)))

ggplot(diff_df, aes(diff, fill = mesh, color=col)) + 
    stat_density(position="identity",geom="line", size=1.5) + scale_color_manual(values=colors[c(1,3,4,2)])
##[which(diff_df$diff>(-400)),]
#geom_density(alpha = 0.2)
#+ scale_fill_manual(values=colors)
#geom_histogram(alpha= 0.5, binwidth = 10, position = 'identity')

year<-2089
ba.file <- paste(ipm.10y.path,"/ba_",(year-9),"_",SIM_NUMBER, ".rdata",sep="")
load(ba.file)
final_ba_10y <- apply(ba,1,sum)#[c(9,10,25,30,77,79,81,90,93,94,96)]

index<-1
len_mesh <- c("m2000", "m3000", "m4000", "m5000")
diff_df <- data.frame(diff=double(), mesh=integer())
for(scn.name in scn.names){
  ba.file <- paste0("./", scn.name, "/ba_",scn.name, "_", iyear, "_", "run_",irun, ".rdata")
  load(ba.file)
  final_ba_1y <- apply(ba,1,sum)#[c(9,10,25,30,77,79,81,90,93,94,96)]
  #ba_diff <- data.frame(diff=(final_ba_1y/final_ba_10y*100-100))#(final_ba_10y - final_ba_1y)
  ba_diff <- data.frame(diff=(final_ba_1y-final_ba_10y))
  rmse_vector[index]<- rmse(final_ba_10y, final_ba_1y)
  ba_diff$mesh <- len_mesh[index]
  diff_df <- rbind(diff_df, ba_diff)
  #hist(ba_diff, breaks=35, main=paste0("ba diff at year ", year, scn.name))
  index<-index+1
}

##mean, median, variance, +-5, +-10, +-20
summary_df<- data.frame(mean=rep(0,4), median= rep(0,4), variance= rep(0,4),
                        a5 = rep(0,4), a10 = rep(0,4), a20= rep(0,4))
index <- 1
for (mesh in len_mesh){
  diff <- diff_df[which(diff_df$mesh==mesh),1]
  summary_df$mean[index]<- round(mean(diff),2)
  summary_df$median[index]<- round(median(diff),2)
  summary_df$variance[index]<- round(var(diff),2)
  summary_df$a5[index]<- length(which(diff<5 & diff>-5))
  summary_df$a10[index]<- length(which(diff<10 & diff>-10))
  summary_df$a20[index]<- length(which(diff<20 & diff>-20))
  index <- index+1
}

write.table(summary_df, file="./summary_basal_area_no_percent.csv", dec=".", sep=";")

ggplot(diff_df, aes(diff, fill = mesh)) + scale_fill_manual(values=colors)+ 
  geom_density(alpha = 0.2)







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