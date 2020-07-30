

spp_comp$IPM.spp <- apply(ba,1, function(x) {ifelse(sum(x)==0,0, which.max(x))})
spp_comp$IPM.spp <- sapply(spp_comp$IPM.spp, function(x) {ifelse(x==0,14,IPM.index.Medfire.spp[x])} )

length(which(spp_comp$spp==spp_comp$IPM.spp))

spp_comp$IPM.spp2 <- apply(ba,1, function(x) {a<-x;a[which.max(x)]=0;ifelse(sum(a)==0,0, which.max(a)) })
spp_comp$IPM.spp2 <- sapply(spp_comp$IPM.spp2, function(x) {ifelse(x==0,14,IPM.index.Medfire.spp[x])} )

length(which((spp_comp$spp==spp_comp$IPM.spp2)&(spp_comp$spp!=14)))

length(which((spp_comp$spp==14) & (spp_comp$IPM.spp== 14)))
length(which(spp_comp$IPM.spp== 14))