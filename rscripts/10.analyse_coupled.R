rm(list=ls())
gc()
setwd("C:/Users/uriso/Desktop")
library(ggplot2)

n.intervals.mesh <- 4000
study.area <- "BCN"
load(paste0("ONEKM/mdl_interface/IPM_",study.area, "_map.rdata"))
orig.adult.trees.file <- paste("ONEKM/IPM/initial_variables/trees_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
orig.ba.file <- paste("ONEKM/IPM/initial_variables/ba_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
orig.saplings.file <- paste("ONEKM/IPM/initial_variables/saplings_", study.area,"_", n.intervals.mesh, ".rdata",sep="")
orig.plots.age.file <- paste("ONEKM/IPM/initial_variables/orig_plots_age", study.area, ".rdata",sep="")
species <- c("phalepensis", "pnigra", "ppinea", "psylvestris", "ppinaster", "puncinata",
             "aalba", "qilex", "qsuber", "qfaginea", "qhumilis", "fsylvatica", "other")
colors <- c("#FFFF99", "#CCFF66", "#99FF33", "#66CC00", "#339900")
load(orig.ba.file); load(orig.saplings.file)
load("ONEKM/Medfire/inputlyrs/rdata/land.rdata")
ini_ba <- ba
ini_sapl<- saplings
ini_land <- land

fire.regeneration <- c(T,T,T,F,T,F,T,F,F,F,T,T,T,T,T,T)
Medfire.index.IPM.spp<-c(5,6,8,9,7,10,2,12,15,11,11,3,1) #quercus humilis and faginea are classified as the same for IPM
IPM.index.Medfire.spp <-c(13,13,12,13,1,2,5,3,4,6,11,8,13,13,9,13)
scenarios <- c("only_IPM_scn_1_run_2", "coupled_scn_1_run_5", "coupled_scn_1_run_6", "coupled_scn_1_run_8")


##mean ba evolution
iyears <- seq(2000,2089,1)
ba_df <- data.frame(matrix(ncol = length(iyears)+1, nrow = length(scenarios)))
colnames(ba_df) <- c("IFN3",iyears) 
rownames(ba_df) <- scenarios
ba_df$IFN3<- mean(apply(ini_ba,1,sum))
index <- 1
for (scenario in scenarios){
  for (iyear in 2000:2089){
    ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
    load(ba.file)
    tot_ba <- apply(ba,1,sum)
    ba_df[index, as.character(iyear)] <- mean(tot_ba)
  }
  index <- index+1
}
plot(seq(2000,2090,1), ba_df[1,], ylim= c(0,25), type = "l", lwd=2, xlab="years", ylab="total basal area")
lines(seq(2000,2090,1),colMeans(ba_df[2:length(scenarios),]), lty= 3, lwd=2)
legend("topleft", legend=c("IPM + MEDFIRE", "IPM only"), lty=c(3,1), cex=0.8)

sd_df <- matrix(ncol = length(iyears)+1, nrow =2)
sd_df[1,]<- 0
sd_df[2,]<- apply(ba_df[2:length(scenarios),],2,sd)

ba_mean_df <- data.frame(years=2000:2090)
ba_mean_df$IPM <- as.double(ba_df[1,])
ba_mean_df$coupled <- colMeans(ba_df[2:length(scenarios),])

ggplot(ba_mean_df, aes(years)) +
  geom_line(aes(y = IPM),color = "purple", size=1.5) + 
  geom_line(aes(y = coupled),  colour = "blue", size=1.5) +
  geom_ribbon(aes(ymin = coupled - sd_df[2,],
                  ymax = coupled + sd_df[2,]), alpha = 0.4, fill="blue") +
  ylim(0,21) + ylab(expression(basal ~ area ~ (m^{2} ~ ha^{-1}))) +
  scale_colour_manual(values=c(colors[4], "yellow"), labels=c("IPM", "COUPLED"))


##num species hist
num_spp_counts <- data.frame(matrix(ncol = 5, nrow = 3))
colnames(num_spp_counts) <- c(0:3,">3") 
rownames(num_spp_counts) <- c("inital data", scenarios)
#num_spp_counts <- integer(7)
for (s in 0:2){
  if(s==0){
    #num_spp <- apply(ini_ba,1,function(x){length(which(x>0))})
    iyear <-2000
    scenario<-scenarios[1]
    ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
    load(ba.file)
    num_spp <- apply(ba,1,function(x){length(which(x>0))})
  } 
  else{
    iyear <-2089
    scenario<-scenarios[s]
    ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
    load(ba.file)
    num_spp <- apply(ba,1,function(x){length(which(x>0))})
  }

  for (i in 0:4){
    if (i!=4){
      num_spp_counts[s+1,i+1]<- length(which(num_spp==i))
    }
    else{
      num_spp_counts[s+1, i+1]<- length(which(num_spp>=i))
    }
  }
}
barplot(t(num_spp_counts), names.arg=c("2000", "2090 - coupled", "2090 - IPM"), col=colors, ylab="number of plots" )
#legend.text=T


num_spp_fractions <- matrix(ncol = 5, nrow = 10)
iyear <-2000
scenario<-scenarios[1]
ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
load(ba.file)
saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear, "_run_1.rdata",sep="")
load(saplings.file)
ini_num_spp <- apply(ba,1,function(x){length(which(x>0))})
ini_num_spp[which(ini_num_spp==0)] <- apply(saplings[which(ini_num_spp==0),],1,function(x){length(which(x>0))})
for (i in 0:4){
    if (i!=4){
      index_num_spp <- which(ini_num_spp==i)
	    length_num_spp <- length(which(ini_num_spp==i))
    } else{
      index_num_spp <- which(ini_num_spp>=i)
	    length_num_spp <- length(which(ini_num_spp>=i))
    }
	for (s in 1:2){
		iyear <-2089
		scenario<-scenarios[s]
		ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
		load(ba.file)
		saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear, "_run_1.rdata",sep="")
		load(saplings.file)
		num_spp <- apply(ba,1,function(x){length(which(x>0))})
		num_spp[which(num_spp==0)] <- apply(saplings[which(num_spp==0),],1,function(x){length(which(x>0))})
		for (j in 0:4){
		    if (j!=4){
			  fraction_num_spp <- length(which(num_spp[index_num_spp]==j))/length_num_spp
		    }
		    else{
			  fraction_num_spp <- length(which(num_spp[index_num_spp]>=j))/length_num_spp
		    }
		    num_spp_fractions[(i+1)+(s-1)*5, j+1] <- fraction_num_spp
  		}
	}
}


interleave <- function(v1,v2)
{
  ord1 <- 2*(1:length(v1))-1
  ord2 <- 2*(1:length(v2))
  c(v1,v2)[order(c(ord1,ord2))]
}
dat<-t(num_spp_fractions[1:5,])
colnames(dat) <- c("0 IPM","1 IPM","2 IPM","3 IPM",">3 IPM")

dat1<-t(num_spp_fractions[6:10,])
colnames(dat1) <- c("0 COUPLED","1 COUPLED","2 COUPLED","3 COUPLED",">3 COUPLED")

dat2 = cbind(dat, dat1)[ , interleave(colnames(dat), colnames(dat1))]

par(mfrow=c(1,5))
lapply(seq(1,ncol(dat2),2), function(x) {
  if(x==1) {
    par(mar=c(5,3,2,0))
    barplot(dat2[,c(x,x+1)], col=colors, xaxt= "n")
  } else {
    par(mar=c(5,2,2,0))
    barplot(dat2[,c(x,x+1)], col=colors, yaxt="n", xaxt= "n")
  }
})


##pinus halapensis %
iyear <-2000
scenario<-scenarios[1]
ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
load(ba.file)
ini_ba <- ba

dom_counts <- data.frame(matrix(ncol = 4, nrow = 3))
dominance_tag<-function(x){
  if (sum(x)==0){
    return(0)
  }
  dom <- x[5]/sum(x)
  if (dom>0.6){
    return(3)
  }
  else if (dom >0.3){
    return(2)
  }
  else if (dom>0){
    return(1)
  }
  else{
    return(0)
  }
}

#dom_dist <- apply(ini_ba[apply(ini_ba,1,sum)!=0,],1,dominance_tag)
for (s in 0:2){
  if(s==0)
    dom_dist <- apply(ini_ba[apply(ini_ba,1,sum)!=0,],1,dominance_tag)
  else{
    iyear <-2089
    scenario<-scenarios[s]
    ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
    load(ba.file)
    saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear, "_run_1.rdata",sep="")
    load(saplings.file)
    dom_dist <- apply(ba[apply(ini_ba,1,sum)!=0,],1,dominance_tag)
    sap_forest <- saplings[apply(ini_ba,1,sum)!=0,]
    dom_dist[which(dom_dist==0)] <- apply(sap_forest[which(dom_dist==0),],1,function(x){ifelse(x[5]>0,1,0)})
    
  }
  
  for (i in 0:3){
    dom_counts[s+1,i+1]<- length(which(dom_dist==i))
  }}
par(mfrow=c(1,1))
barplot(t(dom_counts),  names.arg=c("2000", "2090 - coupled", "2090 - IPM"), col=colors[1:4], ylab="number of plots" )
#legend.text=c("0", "<30%", "30-60%", ">60%"),

dom_matrix <- (matrix(ncol=4, nrow=8))
iyear <-2000
scenario<-scenarios[1]
ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
load(ba.file)
ini_ba <- ba
saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear, "_run_1.rdata",sep="")
load(saplings.file)
ini_dom_dist <- apply(ini_ba[apply(ini_ba,1,sum)!=0,],1,dominance_tag) ##at 2000 there are no forests with only saplings, so ok
sap_forest <- saplings[apply(ini_ba,1,sum)!=0,]
ini_dom_dist[which(ini_dom_dist==0)] <- apply(sap_forest[which(ini_dom_dist==0),],1,function(x){ifelse(x[5]>0,1,0)})

for (i in 0:3){
    
	index_num_spp <- which(ini_dom_dist==i)
	length_num_spp <- length(which(ini_dom_dist==i))

	for (s in 1:2){
		iyear <-2089
		scenario<-scenarios[s]
		ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
		load(ba.file)
		saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear, "_run_1.rdata",sep="")
		load(saplings.file)
		num_spp <- apply(ba[apply(ini_ba,1,sum)!=0,],1, dominance_tag)
		sap_forest <- saplings[apply(ini_ba,1,sum)!=0,]
		num_spp[which(num_spp==0)] <- apply(sap_forest[which(num_spp==0),],1,function(x){ifelse(x[5]>0,1,0)})
		for (j in 0:3){
			fraction_num_spp <- length(which(num_spp[index_num_spp]==j))/length_num_spp
		  dom_matrix[(i+1)+(s-1)*4, j+1] <- fraction_num_spp
  		}
	}
}
dat<-t(dom_matrix[5:8,])
dat1<-t(dom_matrix[1:4,])
dat2 = cbind(dat, dat1)[ , order(c(2*(1:4)-1,2*(1:4)))]

par(mfrow=c(1,4))
lapply(seq(1,ncol(dat2),2), function(x) {
  if(x==1) {
    par(mar=c(5,3,2,0))
    barplot(dat2[,c(x,x+1)], col=colors[-5], xaxt= "n")
  } else {
    par(mar=c(5,2,2,0))
    barplot(dat2[,c(x,x+1)], col=colors[-5], yaxt="n", xaxt= "n")
  }
})



##land evolution

dom_spp_dist <- matrix(ncol = 14, nrow = 4)
dom_spp_dist_sd <- matrix(ncol = 14, nrow = 4)

##land 2010 without shrub
load("./ONEKM/Medfire/inputlyrs/rdata/land.rdata")
ini_land <- land
ini_shrub_bcn_cell_id <- land$cell.id[land$spp==14 & land$cell.id %in% map$Medfire.id]
iyear <- 2010
scenario <- scenarios[2]
land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear, "_run_1.rdata",sep="")
load(land.file)
ini_spp_dist <- hist(land$spp[land$cell.id %in% map$Medfire.id & !(land$cell.id %in% ini_shrub_bcn_cell_id)],  breaks= seq(0,14,1), plot=F)
ini_spp_dist_counts <- ini_spp_dist$counts
dom_spp_dist[1,]<- ini_spp_dist_counts
dom_spp_dist_sd[1,]<- 0

##Coupled
iyear <- 2089

COUPLED_runs <- scenarios[2:length(scenarios)]

spp_dist_counts <- matrix(nrow=length(COUPLED_runs), ncol=14)
index<-1
for (scenario in COUPLED_runs){
  land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear, "_run_1.rdata",sep="")
  load(land.file)
  spp_dist <- hist(land$spp[land$cell.id %in% map$Medfire.id & !(land$cell.id %in% ini_shrub_bcn_cell_id)],breaks= seq(0,14,1), plot=F)
  spp_dist_counts[index,] <- spp_dist$counts
  index <- index+1
}

if (length(COUPLED_runs)>1){
  dom_spp_dist[4,]<- apply(spp_dist_counts,2,mean)
  dom_spp_dist_sd[4,]<- apply(spp_dist_counts,2,sd)
  
} else{
  dom_spp_dist[4,]<- spp_dist_counts
  dom_spp_dist_sd[4,]<- 0
}

##IPM
iyear <- 2089
scenario <- scenarios[1]
ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
load(ba.file)
IPM.to.MEDFIRE.most.abundant<- function(x){
  if (sum(x)==0){
    return(14)
    }
  else{
    spp <- IPM.index.Medfire.spp[which.max(x)] 
    return(spp)}
}
ba_dom_Medfire <- apply(ba,1,IPM.to.MEDFIRE.most.abundant)
ba_dom_Medfire <- ba_dom_Medfire[!(map$Medfire.id %in% ini_shrub_bcn_cell_id)]
spp_dist <- hist(ba_dom_Medfire,breaks= seq(0,14,1), plot=F)
spp_dist_counts <- spp_dist$counts
dom_spp_dist[2,]<- spp_dist_counts
dom_spp_dist_sd[2,]<- 0



##Medfire
iyear <- 2089

MEDFIRE_runs <- c("only_Medfire_scn_1_run_1", "only_Medfire_scn_1_run_2", "only_Medfire_scn_1_run_3",
                  "only_Medfire_scn_1_run_4", "only_Medfire_scn_1_run_5")

spp_dist_counts <- matrix(nrow=length(MEDFIRE_runs), ncol=14)
index<-1
for (scenario in MEDFIRE_runs){
  land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear, "_run_1.rdata",sep="")
  load(land.file)
  spp_dist <- hist(land$spp[land$cell.id %in% map$Medfire.id & !(land$cell.id %in% ini_shrub_bcn_cell_id)],breaks= seq(0,14,1), plot=F)
  spp_dist_counts[index,] <- spp_dist$counts
  index <- index+1
}

if (length(MEDFIRE_runs)>1){
  dom_spp_dist[3,]<- apply(spp_dist_counts,2,mean)
  dom_spp_dist_sd[3,]<- apply(spp_dist_counts,2,sd)
  
} else{
  dom_spp_dist[3,]<- spp_dist_counts
  dom_spp_dist_sd[3,]<- 0
}


# scenario <- "only_Medfire_scn_1_run_1"
# land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear, "_run_1.rdata",sep="")
# load(land.file)
# spp_dist <- hist(land$spp[land$cell.id %in% map$Medfire.id & !(land$cell.id %in% ini_shrub_bcn_cell_id)],breaks= seq(0,14,1), plot=F)
# spp_dist_counts <- spp_dist$counts
# dom_spp_dist[4,]<- spp_dist_counts


dom_spp_dist_barplot <- barplot(dom_spp_dist[,-c(7,10)],  names.arg=c(species[-c(7,10)],"shrub"), cex.names= 0.8, ylim=c(0,max(dom_spp_dist[-c(7,10)])*1.1),
 legend.text=c("2010", "2090 - IPM-1y", "2090 - MEDFIRE-1km", "2090 - COUPLED"),las=2, ylab="number of plots" , col=c("pink","purple", "cyan", "blue"), beside=T)
arrows(dom_spp_dist_barplot, dom_spp_dist[,-c(7,10)]-dom_spp_dist_sd[,-c(7,10)],
       dom_spp_dist_barplot, dom_spp_dist[,-c(7,10)]+dom_spp_dist_sd[,-c(7,10)], angle=90, code=3, length=0.03 )
#legend.text=c("IFN3", "2090 - coupled", "2090 - IPM"),
#barplot(dom_spp_dist, beside=T)
###############
##post-fire analysis
MEDFIRE_runs <- c("only_Medfire_scn_1_run_1", "only_Medfire_scn_1_run_2", "only_Medfire_scn_1_run_3",
                  "only_Medfire_scn_1_run_4", "only_Medfire_scn_1_run_5")
post.fire.mat <- list()
for (i in MEDFIRE_runs)
  post.fire.mat[[i]] <- matrix(ncol = 14, nrow = 14)


for (scenario in MEDFIRE_runs){
  post.fire.succ<- list()
  for (i in 1:14){
    post.fire.succ[[i]]<-integer(0)
  }
  for (iyear in 2010:2089) {
    land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear-1, "_run_1.rdata",sep="")
    load(land.file)
    prev_land <- land
    land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear, "_run_1.rdata",sep="")
    load(land.file)
    burnt.cells <- land$cell.id[(land$tburnt - prev_land$tburnt) == 1 & !is.na(land$tburnt - prev_land$tburnt) & land$cell.id %in% map$Medfire.id]
  
    for (i in 1:14){
      spp_burnt <- which(prev_land$spp==i & land$cell.id %in% burnt.cells)
      post.fire.succ[[i]] <- c(post.fire.succ[[i]], land$spp[spp_burnt])
    }  
  }
  
  
  for (i in 1:14){
    if (length(post.fire.succ[[i]])!=0){
      for (j in 1:14){
        post.fire.mat[[scenario]][i,j]<- length(which(post.fire.succ[[i]]==j))/length(post.fire.succ[[i]])*100
      }
    }
    else{
      post.fire.mat[[scenario]][i,]<-0
    }
  }
  print(paste0("senario: ", scenario," done" ))
}

Y <- do.call(cbind, post.fire.mat)
Y <- array(Y, dim=c(dim(post.fire.mat[[1]]), length(post.fire.mat)))
post.fire.mat <- apply(Y, c(1, 2), mean, na.rm = TRUE)
rownames(post.fire.mat)<- c(species, "shrub")
colnames(post.fire.mat)<- c(species, "shrub")
post.fire.mat[post.fire.mat==0]<-NA
library('plot.matrix')
plot(post.fire.mat[c(-7,-10), c(-7,-10)],fmt.cell='%.2f', las=2,xlab="",
     ylab="" , cex.axis=0.7,cex=0.7, breaks=c(0.00,2,25, 50, 75, 100),
     col=colors, na.print=FALSE)
plot(post.fire.mat[c(-7,-10), c(-7,-10)],fmt.cell='%.2f', axis.col=3, las=2,xlab="",
     ylab="" , cex.axis=0.7,cex=0.7, breaks=c(0.00,2,25, 50, 75, 100),
     col=colors, na.print=FALSE, main="")
#
# post.fire.mat.coupled <- post.fire.mat


##### coupled
MEDFIRE_runs <- scenarios[2:length(scenarios)]
post.fire.mat <- list()
for (i in MEDFIRE_runs)
  post.fire.mat[[i]] <- matrix(ncol = 14, nrow = 14)


for (scenario in MEDFIRE_runs){
  post.fire.succ<- list()
  for (i in 1:14){
    post.fire.succ[[i]]<-integer(0)
  }
  for (iyear in 2010:2089) {
    land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear-1, "_run_1.rdata",sep="")
    load(land.file)
    prev_land <- land
    land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear, "_run_1.rdata",sep="")
    load(land.file)
    burnt.cells <- land$cell.id[(land$tburnt - prev_land$tburnt) == 1 & !is.na(land$tburnt - prev_land$tburnt) & land$cell.id %in% map$Medfire.id]
    
    for (i in 1:14){
      spp_burnt <- which(prev_land$spp==i & land$cell.id %in% burnt.cells)
      post.fire.succ[[i]] <- c(post.fire.succ[[i]], land$spp[spp_burnt])
    }  
  }
  
  
  for (i in 1:14){
    if (length(post.fire.succ[[i]])!=0){
      for (j in 1:14){
        post.fire.mat[[scenario]][i,j]<- length(which(post.fire.succ[[i]]==j))/length(post.fire.succ[[i]])*100
      }
    }
    else{
      post.fire.mat[[scenario]][i,]<-0
    }
  }
  print(paste0("senario: ", scenario," done" ))
}

Y <- do.call(cbind, post.fire.mat)
Y <- array(Y, dim=c(dim(post.fire.mat[[1]]), length(post.fire.mat)))
post.fire.mat <- apply(Y, c(1, 2), mean, na.rm = TRUE)
rownames(post.fire.mat)<- c(species, "shrub")
colnames(post.fire.mat)<- c(species, "shrub")
post.fire.mat[post.fire.mat==0]<-NA
library('plot.matrix')
plot(post.fire.mat[c(-7,-10), c(-7,-10)],fmt.cell='%.2f', axis.col=3, las=2,xlab="",
     ylab="" , cex.axis=0.7,cex=0.7, breaks=c(0.00,0.02,0.25, 0.50,0.75,1),
     col=colors, na.print=FALSE, main="")
# post.fire.mat.coupled <- post.fire.mat
# scenario <- scenarios[1]
# post.fire.succ<- list()
# for (i in 1:14){
#   post.fire.succ[[i]]<-integer(0)
# }
# for (iyear in 2010:2089) {
#   land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear-1, "_run_1.rdata",sep="")
#   load(land.file)
#   prev_land <- land
#   land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear, "_run_1.rdata",sep="")
#   load(land.file)
#   burnt.cells <- land$cell.id[(land$tburnt - prev_land$tburnt) == 1 & !is.na(land$tburnt - prev_land$tburnt) & land$cell.id %in% map$Medfire.id]
# 
#   for (i in 1:14){
#     spp_burnt <- which(prev_land$spp==i & land$cell.id %in% burnt.cells)
#     post.fire.succ[[i]] <- c(post.fire.succ[[i]], land$spp[spp_burnt])
#   }
# }
# 
# post.fire.mat <- matrix(ncol = 14, nrow = 14)
# for (i in 1:14){
#   if (length(post.fire.succ[[i]])!=0){
#     for (j in 1:14){
#       post.fire.mat[i,j]<- length(which(post.fire.succ[[i]]==j))/length(post.fire.succ[[i]])
#     }
#   }
#   else{
#     post.fire.mat[i,]<-0
#   }
# }
plot(post.fire.mat[c(-7,-10), c(-7-10)],fmt.cell='%.3f', cex=0.7, breaks=c(0.01,0.2,0.4,0.6,0.8,1), col=colors)


barplot(t(post.fire.mat[-c(7,10),]), names.arg=c(species[-c(7,10)],"shrub"), col=rainbow(14), ylab="fraction", las=2 )
# legend.text=,
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend = c(species[-c(7,10)], "shrub"), pch=15, pt.cex=1.5, cex=0.9, bty='n',
       col = rainbow(12))




##plot ba hist
years <- c("2000-IFN3", "2090-coupled", "2090-only-IPM")
ba_dist_df <- data.frame(ba=double(), year_model=integer())
aux <- data.frame(ba=ini_tot_ba)
aux$year_model <- years[1]
ba_dist_df <- rbind(ba_dist_df, aux)
aux <- data.frame(ba=tot_ba)
aux$year_model <- years[2]
ba_dist_df <- rbind(ba_dist_df, aux)
# for (y in years[-1]){
#   aux <- data.frame(ba=)
# }
ggplot(ba_dist_df, aes(ba, fill = year_model)) + 
  geom_freqpoly(alpha = 0.2, binwidth=3)


ggplot(ba_dist_df, aes(ba, fill = year_model)) +
  stat_density(position="identity",geom="line")

a <- density(ini_tot_ba, adjust=0.6)
to_ignore <- which(a$x<0)
plot(a$x[-to_ignore], a$y[-to_ignore]*4914, type="l")
b <- density(tot_ba, adjust=0.6)
to_ignore <- which(b$x<0)
lines(b$x[-to_ignore], b$y[-to_ignore]*4914, col="green")

land.file <- paste("remote_output/", scenario,"/land_", scenario,"_", iyear, "_run_1.rdata",sep="")
load(land.file)
IPM.raster <- raster()
extent(IPM.raster)<-extent(c(min(map$X_UTM)*0.99, max(map$X_UTM)*1.01, min(map$Y_UTM)*0.995, max(map$Y_UTM)*1.005)) #c(250000, 540000, 4480000, 4760000)
res(IPM.raster)<-c(1000,1000) #1km^2
crs(IPM.raster) <- CRS("+init=epsg:23030")
IPM.raster <-rasterize(map[,c('X_UTM', 'Y_UTM')], IPM.raster, apply(ba,1,sum))

###PLOT RASTERS
load("./ONEKM/Medfire/inputlyrs/rdata/coordinates.rdata")
library(rgdal)
capitals <- rgdal::readOGR(dsn="./remote_output/CapitalsComarca",layer="CapitalsComarca")

my.id <- match(coord$cell.id, map$Medfire.id)
my.id <- my.id[!is.na(my.id)]


iyear <- 2000
ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
load(ba.file)
ini_ba <- ba


iyear <- 2089
scenario<-scenarios[1]
ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
load(ba.file)
saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear, "_run_1.rdata",sep="")
load(saplings.file)

mean_ba <- ba
mean_ba[,]<-0
mean_sap <- saplings
mean_sap[,]<-0
for (scenario in scenarios[2:length(scenarios)]){
  iyear <- 2089
  ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
  load(ba.file)
  mean_ba <- ba + mean_ba
  saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear, "_run_1.rdata",sep="")
  load(saplings.file)
  mean_sap <- saplings + mean_sap
}
mean_ba <- mean_ba/(length(scenarios)-1)
mean_sap <- mean_sap/(length(scenarios)-1)
ba<- mean_ba
saplings <- mean_sap

IPM.raster <- raster()
extent(IPM.raster)<-c(340000, 500000, 4550000, 4700000)
res(IPM.raster)<-c(1000,1000) #1km^2
crs(IPM.raster) <- CRS("+init=epsg:25831")
##biomass
##IPM.raster <-rasterize(coord[coord$cell.id %in% map$Medfire.id,c(2,3)], IPM.raster, apply(ba[my.id,],1,sum))
# values <- apply(ba,1,sum)
# values[which(values==0)] <- apply(saplings[which(values==0),],1,function(x){ifelse(length(which(x>0))>0,0.1,0)})
# IPM.raster <-rasterize(coord[coord$cell.id %in% map$Medfire.id,c(2,3)], IPM.raster, values[my.id])

##spp count
# values <- apply(ba,1,function(x){length(which(x>0))})
# values[which(values==0)] <- apply(saplings[which(values==0),],1,function(x){length(which(x>0))})

mean_values <- double(nrow(ba))
for (scenario in scenarios[2:length(scenarios)]){
  iyear <- 2089
  ba.file <- paste("remote_output/", scenario,"/ba_", scenario,"_", iyear, "_run_1.rdata",sep="")
  load(ba.file)
  saplings.file <- paste("remote_output/", scenario,"/saplings_", scenario,"_", iyear, "_run_1.rdata",sep="")
  load(saplings.file)
  values <- apply(ba,1,function(x){length(which(x>0))})
  values[which(values==0)] <- apply(saplings[which(values==0),],1,function(x){length(which(x>0))})
  mean_values <- values + mean_values
}
mean_values <- mean_values/(length(scenarios)-1)
values<-mean_values
IPM.raster <-rasterize(coord[coord$cell.id %in% map$Medfire.id,c(2,3)], IPM.raster, values[my.id])



##BCN province maps
#b <- which(coord$y>4550000 & coord$y<4700000 & coord$x>350000 & coord$x<500000)
BCN.raster <- raster()
extent(BCN.raster)<-c(340000, 500000, 4550000, 4700000)
res(BCN.raster)<-c(1000,1000) #1km^2
crs(BCN.raster) <- CRS("+init=epsg:25831")
BCN.raster <-rasterize(coordinates(coord[coord$cell.id %in% cell.id.interface.bcn$Medfire.id,c(2,3)]), BCN.raster, -30) #land$tburnt[land$cell.id %in% cell.id.interface.bcn$Medfire.id]

to.fill.index <- which(!is.na(IPM.raster[]))
BCN.raster[to.fill.index] <- IPM.raster[to.fill.index]
##Note that when giving breakpoints, the default for R is that the histogram 
##cells are right-closed (left open) intervals of the form (a,b].
##You can change this with the right=FALSE option, 
##which would change the intervals to be of the form [a,b). 


#cuts<-c(-31,-0.01,0.01,8,16,24,1000)
cuts<-c(-31,-0.01,0.01,1.5,2.5,3.5,16)
#cuts<-c(-31,-29,-2.1,-1.1, -0.1, 0.1, 1.9, 2.9, 16, 50)
pal <- colorRampPalette(c("orange", "yellow"))
colors <- c("#FFFF99", "#CCFF66", "#99FF33", "#66CC00", "#339900")
#plot(BCN.raster, breaks=cuts, col=c("gray", pal(3) , colors[-1], "red"))
plot(BCN.raster, breaks=cuts, col=c("gray",colors))
plot(capitals.bcn, pch=17, bg="transparent", add=TRUE)
posit <- c(4,4,4,3,4,1,4,4,4,4,4,4)
text(capitals.bcn, labels = 1:13, pos = 4, offset = 0.2, cex=0.7)
scalebar(40000, xy=c(440000, 4560000), type='bar', divs=4, label=c(0, 20, 40), below="Km", cex=0.5)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# legend("topleft", legend = c( "urban/crops", "shrub", "0-8", "8-16", "16-24",">24","capital comarca" ), pch=c(17, rep(15,6)), pt.cex=c(2, rep(3,6)), cex=1.5, bty='n',
#        col = c("black","gray",colors))
legend("topleft", legend = c( "urban/crops", "shrub", "0-8", "8-16", "16-24",">24","capital comarca" ), pch=c(rep(15,6),17), pt.cex=c(1, rep(1,6)), cex=0.7, bty='n',
       col = c("gray",colors,"black"), ncol =4)
legend("topleft", legend = c("urban/crops","shrub","0-1.5","1.5-2.5","2.5-3.5",">3.5", "capital comarca" ), pch=c(rep(15,6),17), pt.cex=c(1, rep(1,6)), cex=0.7, bty='n',
       col = c("gray",colors,"black"), ncol =4)
legend("topleft", legend = c("0", "<30%", "30-60%", ">60%"), pch=15, pt.cex=3, cex=1.5, bty='n',
       col = c(colors[-5]))
legend("topleft", legend = c("shrub","<-2", "-2", "-1", "0", "1", "2", ">2"), pch=15, pt.cex=2, cex=1, bty='n',
       col = c("red", pal(3) , colors[-1]))
legend("topleft", legend = c("IPM-1y","COUPLED"), pch=15, pt.cex=3, cex=1.5, bty='n',
       col = c("purple", "blue"))
