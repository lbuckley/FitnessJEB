library(ggplot2)
library(cowplot)
library(reshape2)

years= 1950:2099
#absorptivity
aseq= seq(0.4,0.7,0.05)

#CHOOSE PROJECTION
proj.k=2 #1: bcc-csm1-1.1.rcp60, 2: ccsm4.1.rcp60, 3: gfdl-cm3.1.rcp60
projs=c("bcc-csm","ccsm4","gfdl")

setwd(wd)
setwd("./Data/")
ldat<-read.csv("ColiasLambdaData.csv")
ldat$year= as.factor(ldat$year)

#plot
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/figures/")
pdf("ButterflyFig.pdf", height = 8, width = 8)
ggplot(data=ldat, aes(x=elev, y=value, color=year))+geom_point()+facet_grid(fitcomp~absorp, scales="free")+
  ylab("fitness component")+xlab("elevation (m)")+geom_smooth(method="loess",se=FALSE, aes(group=year))+theme_bw(base_size = 16)+
  theme(legend.position="bottom")+
  scale_color_viridis(discrete=TRUE) 
dev.off()



