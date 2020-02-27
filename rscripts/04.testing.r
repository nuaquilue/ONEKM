rm(list=ls())
library(tidyverse)
setwd("c:/work/MEDMOD/SpatialModelsR/MEDFIRE")  #N?Laptop
load("inputlyrs/rdata/land.rdata")

ggplot(filter(land, spp<=13, age<=10), aes(x=as.factor(age), y=biom)) +
  geom_violin(width=1.2, color="black", fill="grey") +
  geom_boxplot(width=0.1) +  theme_bw() +   theme(legend.position="none") + 
  facet_wrap(.~as.factor(spp)) +
  stat_summary(fun.y=median, geom="point", size=2, color="black", shape="square") 

ggplot(filter(land, spp<=13, age>10 & age<=20), aes(x=as.factor(age), y=biom)) +
  geom_violin(width=1.2, color="black", fill="grey") +
  geom_boxplot(width=0.1) +  theme_bw() +   theme(legend.position="none") + 
  facet_wrap(.~as.factor(spp)) +
  stat_summary(fun.y=median, geom="point", size=2, color="black", shape="square") 



