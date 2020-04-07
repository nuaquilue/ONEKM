load(paste0(out.path, "/burnt_cells.rdata"))
load("inputlyrs/rdata/pfst.pwind.rdata")
swc <- 1
pigni <- data.frame(cell.id=land$cell.id, p=pigni*pfst.pwind[,ifelse(swc==1,1,2)])
pigni <- filter(pigni, !is.na(p) & p>0)
pfst.pwind$cell.id <- land$cell.id

old.forest.coord <- filter(land, spp<=3 & age>=30) %>% select(cell.id) %>% left_join(coord, by = "cell.id")

igni.id <- sample(pigni$cell.id, 1, replace=F, pigni$p)
fst.sprd.weight <- read.table("inputfiles/SprdRateWeights.txt", header=T)
spp.flammability <- read.table("inputfiles/SppSpreadRate.txt", header=T)
wwind <- fst.sprd.weight[1,fire.spread.type+1]
wslope <- fst.sprd.weight[2,fire.spread.type+1]
wfuel <- fst.sprd.weight[3,fire.spread.type+1]
wflam <- fst.sprd.weight[4,fire.spread.type+1]
waspc <- fst.sprd.weight[5,fire.spread.type+1]
pfst.pwind$cell.id <- land$cell.id
fire.wind <- sample(c(0,315,270), 1, replace=F, p=filter(pfst.pwind,cell.id==igni.id)[3:5])

default.windir <- data.frame(x=c(0,-1,1,290,-290,289,-291,291,-289,-2,2,580,-580),
                             windir=c(-1,270,90,180,0,225,45,135,315,270,90,180,0))
neighs <- nn2(coord[,-1], filter(coord, cell.id %in% fire.front)[,-1], searchtype="priority", k=13)

neigh.id <- data.frame(cell.id=coord$cell.id[neighs$nn.idx],
                       source.id=rep(sort(fire.front), 13),
                       dist=as.numeric(neighs$nn.dists))
neigh.id <- filter(neigh.id, dist<=2000) %>% mutate(x=cell.id-source.id) %>% 
  left_join(default.windir, by="x") %>% select(-x) 

neigh.id <- filter(neigh.id, cell.id %notin% visit.cells)  #before it was burnt.cells

## Filter 'orography' for source and neigbour cells
neigh.orography <- filter(orography, cell.id %in% c(fire.front, neigh.id$cell.id)) %>% select(cell.id, elev, aspect)

## Compute species flammability and aspect factors for spread. rate
flam <- filter(land, cell.id %in% neigh.id$cell.id) %>% select(cell.id, spp) %>%
  left_join(spp.flammability[,c(1,fire.spread.type+1)], by="spp") 
flam$y <- wflam * flam[,ncol(flam)]

aspc <- filter(neigh.orography, cell.id %in% neigh.id$cell.id) %>% select(cell.id, aspect) %>%
  mutate(z=waspc*ifelse(aspect==1, 0.1, ifelse(aspect==3, 0.9, ifelse(aspect==4, 0.4, 0.3))))

sprd.rate <-  left_join(neigh.id, select(neigh.orography, cell.id, elev), by="cell.id") %>%
  left_join(select(neigh.orography, cell.id, elev), by=c("source.id"="cell.id")) %>% 
  mutate(dif.elev = elev.x-elev.y, 
         front.slope = wslope * pmax(pmin(dif.elev/dist,0.5),-0.5)+0.5, 
         front.wind = wwind * (ifelse(abs(windir-fire.wind)>180, 
                                      360-abs(windir-fire.wind), abs(windir-fire.wind)))/180) %>%
  left_join(select(flam, cell.id, y), by="cell.id") %>% 
  left_join(select(aspc, cell.id, z), by="cell.id") %>%
  mutate(sr.noacc=front.slope+front.wind+y+z,
         # sr=round(sr.noacc*fire.strength,3),
         # pb1=(1-exp(-sr.noacc))^rpb,
         # pb2=(1-exp(-sr))^rpb,
         # pb3=1-exp(-sr.noacc^rpb),
         pb4=sr.noacc^rpb,
         pb5=1+rpb*log(sr.noacc)) %>% #select(-front.slope, -front.wind, -y, -z)
  group_by(cell.id) %>% summarize(sr.noacc=max(sr.noacc), pb=max(pb5))

sprd.rate$burning <- runif(nrow(sprd.rate), 0, pb.upper.th) <= sprd.rate$pb & sprd.rate$pb> pb.lower.th #

burnt.cells <- c(burnt.cells, sprd.rate$cell.id[sprd.rate$burning])
visit.cells <- c(visit.cells, sprd.rate$cell.id)
burnt.intens <- c(burnt.intens, sprd.rate$sr.noacc[sprd.rate$burning]>ifelse(swc<4,fire.intens.th,100))
exclude.th <- min(max(sprd.rate$sr.noacc)-0.005, 
                  rnorm(1,mean(sprd.rate$sr.noacc[sprd.rate$burning])-mad(sprd.rate$sr.noacc[sprd.rate$burning])/2,mad(sprd.rate$sr.noacc[sprd.rate$burning])))

fire.front <- sprd.rate$cell.id[sprd.rate$burning & sprd.rate$sr.noacc>=exclude.th]
