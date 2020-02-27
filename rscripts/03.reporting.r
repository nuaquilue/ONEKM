play <- function(){
  rm(list=ls())  
  setwd("c:/work/MEDMOD/SpatialModelsR/MEDFIRE")  #N?Laptop
  setwd("d:/MEDMOD/SpatialModelsR/MEDFIRE")   #CTFC
  
  source("rscripts/03.reporting.r")
  
  scn.name <- "Test03_rcp85_1p"
  plot.abundance(scn.name)
  plot.drought(scn.name)
  plot.cohort(scn.name)
  plot.afforest(scn.name)
  
}

species.name.color <- function(){
  species <- data.frame(spp=1:14, name=c("phalepensis", "pnigra", "ppinea", "psylvestris", "ppinaster",
                                         "puncinata", "aalba", "qilex", "qsuber", "qfaginea", "qhumilis",
                                         "fsylvatica", "other", "shrub"),
                        name.sort=c("aalba","fsylvatica","other","phalepensis","pnigra","ppinaster","ppinea",
                                    "psylvestris", "puncinata", "qfaginea", "qhumilis", "qilex", "qsuber", "shrub"),
                        color=c("chartreuse3", "darkolivegreen1", "darkseagreen", "forestgreen",
                                "olivedrab4", "darkslategrey", "blue4", "gold", "saddlebrown",
                                "sienna2", "palegoldenrod", "red3", "purple3", "grey70"),
                        color.sort=c("royalblue4","red4", "purple3", "chartreuse2", "darkolivegreen1",
                                     "darkslategrey",  "seagreen1", "forestgreen", "grey10",   
                                     "darkorange2", "palegoldenrod", "gold", "saddlebrown", "grey70"),
                        type=c("conifer", "conifer", "conifer", "conifer", "conifer",
                               "conifer", "conifer", "deciduous", "deciduous", "deciduous", "deciduous",
                               "deciduous", "deciduous", "shrub"))
  
  
  save(species, file="rscripts/ins/species.rdata")
}


plot.drought <- function(scn.name){
  
  library(tidyverse)
  load("rscripts/ins/species.rdata")  
  
  ## Existences
  dta.land <- read.table(paste0("outputs/", scn.name, "/Land.txt"), header=T)
  dta.land <- mutate(dta.land, decade= ((dta.land$year-1) %/% 10)*10+2010) %>%
              mutate(decade=ifelse(decade==2100,2090, decade)) %>%
              group_by(decade, spp, run) %>% summarise(area=sum(vol))%>%
              group_by(decade, spp) %>% summarise(area=mean(area)) 
  ## Area and pctg killed by drought
  dta.drought <- read.table(paste0("outputs/", scn.name, "/Drought.txt"), header=T)
  dta.drought <- mutate(dta.drought, decade=((year-1) %/% 10)*10+2010) %>%
                 mutate(decade=ifelse(decade==2100,2090, decade)) %>%
                 group_by(decade, spp, run) %>% summarise(ha=sum(ha))%>%
                 group_by(decade, spp) %>% summarise(ha=mean(ha)) %>% 
                 left_join(dta.land, by=c("decade", "spp")) %>% 
                 left_join(select(species, spp, name), by="spp") %>%
                 mutate(pctg.kill=100*ha/area)
  ## Accumualted area and pctg killed by drought
  dta.accum <- data.frame(decade=2010, spp=1:13)
  dta.accum <- left_join(dta.accum, filter(dta.drought, decade==2010) %>% select(decade, spp, ha), by = c("decade", "spp")) %>%
               left_join(dta.land, by = c("decade", "spp"))
  dta.accum$ha[is.na(dta.accum$ha)] <- 0
  for(d in seq(2020,2090,10)){
    aux <- data.frame(decade=d, spp=1:13)
    aux <- left_join(aux,  filter(dta.drought, decade==d) %>% select(decade, spp, ha), by = c("decade", "spp"))%>%
           left_join(dta.land, by = c("decade", "spp"))
    aux$ha[is.na(aux$ha)] <- 0
    dta.accum <- rbind(dta.accum,aux)
    dta.accum$ha[dta.accum$decade==d] <- dta.accum$ha[dta.accum$decade==d] + dta.accum$ha[dta.accum$decade==(d-10)]
  }
  dta.accum$pctg.kill <- 100*dta.accum$ha/dta.accum$area
  dta.accum$name <- rep(species$name[-14],1)
    
  ## PLOT Area killed by decade (km2) 
  idx <- which(species$name.sort %in% sort(unique(dta.drought$name)))
  p1 <- ggplot(data=dta.drought, aes(x=as.factor(decade), y=ha/100, fill=name)) + geom_bar(stat="identity") +
        scale_fill_manual("Species", values=as.character(species$color.sort[idx])) + 
        ggtitle("Area killed by drought per decade") + ylab("km2") + xlab("period") +
        theme_bw() + theme(axis.text.x = element_text(angle = 90)) 
  ## PLOT Pctg killed over existence by decade (%) 
  p2 <- ggplot(data=dta.drought, aes(x=as.factor(decade), y=pctg.kill, fill=name)) + geom_bar(stat="identity") +
        scale_fill_manual("Species", values=as.character(species$color.sort[idx])) + 
        ggtitle("Pctg of actual area killed by drought per decade") + ylab("%") + xlab("period") +
        theme_bw() + theme(axis.text.x = element_text(angle = 90)) 
  ## PLOT Accumulated area killed by decade (km2) 
  p3 <- ggplot(data=dta.accum, aes(x=as.factor(decade), y=ha/100, fill=name)) + geom_bar(stat="identity") +
        scale_fill_manual("Species", values=as.character(species$color.sort)) + 
        ggtitle("Accumulated area killed by drought per decade") + ylab("km2") + xlab("period") +
        theme_bw() + theme(axis.text.x = element_text(angle = 90)) 
  ## PLOT Pctg killed over existence by decade (%) 
  p4 <- ggplot(data=dta.accum, aes(x=as.factor(decade), y=pctg.kill, fill=name)) + geom_bar(stat="identity") +
        scale_fill_manual("Species", values=as.character(species$color.sort)) + 
        ggtitle("Accumulated pctg of actual area killed by drought per decade") + ylab("%") + xlab("period") +
        theme_bw() + theme(axis.text.x = element_text(angle = 90)) 
  
  ## Save in tiff
  tiff(paste0("rscripts/outs/DroughtKilled_", scn.name, ".tiff"), width=1200, height=800)
  gridExtra::grid.arrange(p1,p2,p3,p4, nrow=2)
  dev.off()
}


plot.cohort <- function(scn.name){
  
  library(tidyverse)
  load("rscripts/ins/species.rdata") 
  
  ## Area of species replaced after drought
  dta.cohort <- read.table(paste0("outputs/", scn.name, "/cohort.txt"), header=T)
  dta.cohort <- group_by(dta.cohort, run, spp.out, spp.in) %>% summarise(ha=sum(ha)) %>%
                group_by(spp.out, spp.in) %>% summarise(ha=mean(ha)) %>% 
                left_join(select(species, spp, name), by=c("spp.out"="spp") ) %>%
                mutate(name.out=name) %>% select(-name) %>%
                left_join(select(species, spp, name), by=c("spp.in"="spp") ) %>%
                mutate(name.in=name) %>% select(-name)
  all.out <- group_by(dta.cohort, spp.out) %>% summarise(tot.ha=sum(ha))
  dta.cohort <- left_join(dta.cohort, all.out, by="spp.out") %>%
                mutate(pctg=100*ha/tot.ha)
    
  ##  PLOT area replaced by new spp / killed spp
  select.color <- filter(species, name %in% unique(dta.cohort$name.in)) %>% select(color.sort)
  p1 <- ggplot(data=dta.cohort, aes(x=name.out, y=ha/100, fill=name.in)) + geom_bar(stat="identity") +
        scale_fill_manual("Species", values=as.character(select.color$color.sort) )+ 
        ggtitle("Species replacement after drought") + ylab("km2") + xlab("species killed") +
        theme_bw() + theme(axis.text.x = element_text(angle = 90)) 
  ## PLOT pctg of area replaced by new spp / killed spp
  p2 <- ggplot(data=dta.cohort, aes(x=name.out, y=pctg, fill=name.in)) + geom_bar(stat="identity") +
        scale_fill_manual("Species", values=as.character(select.color$color.sort) )+ 
        ggtitle("Percentage of species replacement after drought") + ylab("%") + xlab("species killed") +
        theme_bw() + theme(axis.text.x = element_text(angle = 90)) 
 
  ## Save in tiff
  tiff(paste0("rscripts/outs/CohortReplace_", scn.name, ".tiff"), width=1200, height=500)
  gridExtra::grid.arrange(p1,p2, nrow=1)
  dev.off() 
}


plot.abundance <- function(scn.name){
  
  library(tidyverse)
  load("rscripts/ins/species.rdata")  
  
  ## Existences
  dta.land <- read.table(paste0("outputs/", scn.name, "/Land.txt"), header=T)
  dta.land <- mutate(dta.land, year=year+2009)
  dta.ini <- filter(dta.land, year==2010) %>% select(-year)
  names(dta.ini)[3:6] <- paste0(names(dta.ini)[3:6], ".ini") 
  dta.land <- left_join(dta.land, dta.ini, by=c("run","spp")) %>%
              mutate(pctg.area=100*area/area.ini, pctg.vol=100*vol/vol.ini,
                     pctg.volbark=100*volbark/volbark.ini,
                     pctg.carbon=100*carbon/carbon.ini,
                     rel.vol=vol/area, rel.volbark=volbark/area, rel.carbon=carbon/area) %>%
              left_join(select(species, spp, name, type), by="spp")
  
  ## PLOT Abundance area
  p1 <- ggplot(data=filter(dta.land,spp<=13), aes(x=year, y=area/100, colour=name)) +
        geom_smooth(method = "loess", size=1.5) + facet_wrap(.~type) +
        scale_color_manual("Species", values=as.character(species$color.sort)) + 
        ggtitle("Abundance - Area") + ylab("km2") + theme_bw()
  ## PLOT Abundance volume with bark
  p2 <- ggplot(data=filter(dta.land,spp<=13), aes(x=year, y=volbark/10^6, colour=name)) +
        geom_smooth(method = "loess", size=1.5) + facet_wrap(.~type) +
        scale_color_manual("Species",values=as.character(species$color.sort)) + 
        ggtitle("Abundance - Volume with bark") + ylab("m3/10^6") + theme_bw() 
  ## PLOT Abundance carbon
  p3 <- ggplot(data=filter(dta.land,spp<=13), aes(x=year, y=carbon/10^6, colour=name)) +
        geom_smooth( method = "loess", size=1.5) + facet_wrap(.~type) +
        scale_color_manual("Species", values=as.character(species$color.sort)) + 
        ggtitle("Abundance - Carbon") + ylab("t/10^6") + theme_bw()
  ## PLOT Percentage area
  p4 <- ggplot(data=filter(dta.land,spp<=13), aes(x=year, y=pctg.area, colour=name)) +
        geom_smooth(method = "loess", size=1.5) + facet_wrap(.~type) +
        scale_color_manual("Species", values=as.character(species$color.sort)) + 
        ggtitle("Percentage - Area") + ylab("%") + theme_bw()
  ## PLOT Percentage volume with bark
  p5 <- ggplot(data=filter(dta.land,spp<=13), aes(x=year, y=pctg.volbark, colour=name)) +
        geom_smooth(method = "loess", size=1.5) + facet_wrap(.~type) +
        scale_color_manual("Species",values=as.character(species$color.sort)) + 
        ggtitle("Percentage - Volume with bark") + ylab("%") + theme_bw() 
  ## PLOT Percentage carbon
  p6 <- ggplot(data=filter(dta.land,spp<=13), aes(x=year, y=pctg.carbon, colour=name)) +
        geom_smooth( method = "loess", size=1.5) + facet_wrap(.~type) +
        scale_color_manual("Species", values=as.character(species$color.sort)) + 
        ggtitle("Percentage - Carbon") + ylab("%") + theme_bw()
  ## PLOT Relative increment volume with bark
  p7 <- ggplot(data=filter(dta.land,spp<=13), aes(x=year, y=rel.volbark, colour=name)) +
        geom_smooth(method = "loess", size=1.5) + facet_wrap(.~type) +
        scale_color_manual("Species",values=as.character(species$color.sort)) + 
        ggtitle("Relative abundance - Volume with bark") + ylab("m3/ha") + theme_bw() 
  ## PLOT Relative increment carbon
  p8 <- ggplot(data=filter(dta.land,spp<=13), aes(x=year, y=rel.carbon, colour=name)) +
        geom_smooth( method = "loess", size=1.5) + facet_wrap(.~type) +
        scale_color_manual("Species", values=as.character(species$color.sort)) + 
        ggtitle("Relative abundance - Carbon") + ylab("t/ha") + theme_bw()
  
  ## PLOTs for SHRUBS
  p9 <- ggplot(data=filter(dta.land,spp==14), aes(x=year, y=area/100, colour=name)) +
        geom_smooth(method = "loess", size=1.5) + 
        scale_color_manual("Species", values=as.character(species$color.sort)[14]) + 
        ggtitle("Abundance - Area") + ylab("km2") + theme_bw()
  p10 <- ggplot(data=filter(dta.land,spp==14), aes(x=year, y=vol/10^9, colour=name)) +
         geom_smooth(method = "loess", size=1.5) + 
         scale_color_manual("Species",values=as.character(species$color.sort)[14]) + 
         ggtitle("Abundance - Mass") + ylab("t/10^9") + theme_bw() 
  p11 <- ggplot(data=filter(dta.land,spp==14), aes(x=year, y=pctg.area, colour=name)) +
         geom_smooth(method = "loess", size=1.5) +
         scale_color_manual("Species", values=as.character(species$color.sort)[14]) + 
         ggtitle("Percentage - Area") + ylab("%") + theme_bw()
  p12 <- ggplot(data=filter(dta.land,spp==14), aes(x=year, y=pctg.vol, colour=name)) +
         geom_smooth(method = "loess", size=1.5) + 
         scale_color_manual("Species",values=as.character(species$color.sort)[14]) + 
         ggtitle("Percentage - Mass") + ylab("%") + theme_bw() 
  
  
  ## Save in tiff
  tiff(paste0("rscripts/outs/Abundance_", scn.name, ".tiff"), width=1800, height=500)
  gridExtra::grid.arrange(p1,p2,p3, nrow=1)
  dev.off() 
  tiff(paste0("rscripts/outs/AbundancePctg_", scn.name, ".tiff"), width=1800, height=500)
  gridExtra::grid.arrange(p4,p5,p6, nrow=1)
  dev.off() 
  tiff(paste0("rscripts/outs/AbundanceRel_", scn.name, ".tiff"), width=1200, height=500)
  gridExtra::grid.arrange(p7,p8, nrow=1)
  dev.off() 
  tiff(paste0("rscripts/outs/AbundanceShrub_", scn.name, ".tiff"), width=1200, height=500)
  gridExtra::grid.arrange(p9,p10,p11,p12, nrow=2)
  dev.off() 
  
  
}


plot.afforest <- function(scn.name){
  
  library(tidyverse)
  load("rscripts/ins/species.rdata") 
  
  ## Area of species replaced after drought
  dta.afforest <- read.table(paste0("outputs/", scn.name, "/Afforestation.txt"), header=T)
  dta.afforest <- mutate(dta.afforest, decade=((year-1) %/% 10)*10+2010) %>%
                  mutate(decade=ifelse(decade==2100,2090, decade)) %>%
                  group_by(decade, spp, run) %>% summarise(ha=sum(ha))%>%
                  group_by(decade, spp) %>% summarise(ha=mean(ha)) %>% 
                  left_join(select(species, spp, name), by="spp") #%>%
  afforest.decade <- group_by(dta.afforest, decade) %>% summarise(tot.ha=sum(ha))
  dta.afforest <- left_join(dta.afforest, afforest.decade, by="decade") %>%
                  mutate(pctg.ha=100*ha/tot.ha)
  
  ## PLOT Area colonized  by decade (km2) 
  p1 <- ggplot(data=dta.afforest, aes(x=as.factor(decade), y=ha/100, fill=name)) + geom_bar(stat="identity") +
        scale_fill_manual("Species", values=as.character(species$color.sort)) + 
        ggtitle("Area colonized by tree species per decade") + ylab("km2") + xlab("period") +
        theme_bw() + theme(axis.text.x = element_text(angle = 90)) 
  ## PLOT Pctg area colonized over total colonized  by decade (%) 
  p2 <- ggplot(data=dta.afforest, aes(x=as.factor(decade), y=pctg.ha, fill=name)) + geom_bar(stat="identity") +
        scale_fill_manual("Species", values=as.character(species$color.sort)) + 
        ggtitle("Pctg of area colonized by tree species per decade") + ylab("%") + xlab("period") +
        theme_bw() + theme(axis.text.x = element_text(angle = 90)) 
  
  ## Save in tiff
  tiff(paste0("rscripts/outs/Afforestation_", scn.name, ".tiff"), width=1200, height=500)
  gridExtra::grid.arrange(p1,p2, nrow=1)
  dev.off() 
}