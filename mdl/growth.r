######################################################################################
##
######################################################################################

growth <- function(land, clim){
 
  ## Tracking
  cat("Species growth", "\n")

  ## Read coefficients
  growth.coeff <- read.table("inputfiles/GrowthSpp.txt", header=T)
  growth.coeff.shrub <- read.table("inputfiles/GrowthShrub.txt", header=T)
  
  ## Growth of species when SDM in, species when SDM out, and shrub
  aux.spp <- filter(land, spp<=13) %>% select(cell.id, spp, biom) %>% 
             left_join(select(clim, cell.id, sdm, sqi), by = "cell.id") 
  aux.spp.sdmin <- filter(aux.spp, sdm==1) %>% left_join(growth.coeff, by = c("spp", "sqi")) %>% 
                   mutate(increment = x*biom/10 + x2*(biom/10)^2 ) %>%
                   mutate(increment=ifelse(increment<c,c,increment)) %>% select(cell.id, increment)
  aux.spp.sdmout <- filter(aux.spp, sdm==0) %>% left_join(filter(growth.coeff, sqi==1) %>% select(-sqi), by = "spp") %>% 
                    mutate(increment = x*biom/10 + x2*(biom/10)^2 ) %>%
                    mutate(increment=ifelse(increment<c,c,increment)) %>% select(cell.id, increment)
  aux.shrub <- filter(land, spp==14) %>% select(cell.id, spp, biom) %>%
               left_join(select(clim, cell.id, sdm, sqi), by = "cell.id")         
  aux.shrub <- left_join(aux.shrub, growth.coeff.shrub, by = "sqi")  %>%
               mutate(increment = a*log(biom/10000) + b ) %>% 
               mutate(increment = ifelse(increment<=0, b*1000, increment)) %>% select(cell.id, increment)
  
  
  ## Join increment
  all <- rbind(aux.spp.sdmin, aux.spp.sdmout, aux.shrub)
  new.biom <- left_join(land, all, by = "cell.id") %>% mutate(biom=biom+increment*10) %>% select(biom)
  
  return(new.biom$biom)
    
}
  

growth.10y <- function(land, clim){
  
  ## Tracking
  cat("Species growth", "\n")
  
  ## Read coefficients
  growth.coeff <- read.table("inputfiles/GrowthSpp.txt", header=T)
  growth.coeff.shrub <- read.table("inputfiles/GrowthShrub.txt", header=T)
  ba.ini <- read.table("inputfiles/SppBasalArea10y.txt", header=T)
  
  ## Growth of species when SDM in, species when SDM out, and shrub
  aux.spp <- filter(land, spp<=13) %>% select(cell.id, spp, biom, age) %>% 
             left_join(select(clim, cell.id, sdm, sqi), by = "cell.id") 
  aux.spp.sdmin <- filter(aux.spp, sdm==1 & age>10) %>% left_join(growth.coeff, by = c("spp", "sqi")) %>% 
                   mutate(increment = c+x*biom/10 + x2*(biom/10)^2)  %>% 
                   mutate(biom=biom+increment*10) %>% select(cell.id, biom)
  aux.spp.sdmin.10y <- filter(aux.spp, sdm==1 & age<=10) %>% left_join(ba.ini, by = c("spp", "sqi", "age")) %>% 
                       mutate(biom=biom.fix) %>% select(cell.id, biom)
  aux.spp.sdmout <- filter(aux.spp, sdm==0 & age>10) %>% 
                    left_join(filter(growth.coeff, sqi==1) %>% select(-sqi), by = "spp") %>% 
                    mutate(increment = c+x*biom/10 + x2*(biom/10)^2 ) %>%
                    mutate(biom=biom+increment*10) %>% select(cell.id, biom)
  aux.spp.sdmout.10y <- filter(aux.spp, sdm==0 & age<=10) %>% 
                        left_join(filter(ba.ini, sqi==1) %>% select(-sqi), by = c("spp", "age")) %>% 
                        mutate(biom=biom.fix) %>% select(cell.id, biom)
  aux.shrub <- filter(land, spp==14) %>% select(cell.id, spp, biom) %>%
               left_join(select(clim, cell.id, sdm, sqi), by = "cell.id")         
  aux.shrub <- left_join(aux.shrub, growth.coeff.shrub, by = "sqi")  %>%
               mutate(increment = a*log(biom/10000) + b ) %>% 
               mutate(increment = ifelse(increment<=0, b*1000, increment)) %>% 
               mutate(biom=biom+increment*10) %>% select(cell.id, biom)
  
  
  ## Join increment
  all <- rbind(aux.spp.sdmin, aux.spp.sdmin.10y, aux.spp.sdmout, aux.spp.sdmout.10y, aux.shrub)
  new.biom <- select(land, -biom) %>% left_join(all, by = "cell.id")  %>% select(biom)
  
  return(new.biom$biom)
  
}
