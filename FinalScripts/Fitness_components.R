#Explore survival and fecundity data in meta-analyses
library(ggplot2)
library(reshape2)
library(Hmisc)
library(cowplot)
library(viridis)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/meta_data/Scranton_TPC/")
med= read.csv("data_med_R.csv")
temp= read.csv("data_temp_R.csv")
trop= read.csv("data_trop_R.csv")
med$pop="Mediterranean"
temp$pop="temperate"
trop$pop="tropical"

tpc= rbind(med,temp,trop)
tpc$asurv= 1-tpc$amort
tpc$jsurv= 1-tpc$jmort
#to long format
tpc= melt(tpc, id.vars = c("TempK","pop"), measure.vars = c("brate","mat","jmort","amort","fec","plotmat","asurv","jsurv"))
#select fitness components
tpc= tpc[tpc$variable %in% c("brate","mat","asurv"),]
#normalize fecund
tpc[tpc$pop=="Mediterranean"&tpc$variable=="brate","value"]= tpc[tpc$pop=="Mediterranean"&tpc$variable=="brate","value"]/ max(tpc[tpc$pop=="Mediterranean"&tpc$variable=="brate","value"], na.rm=TRUE)
tpc[tpc$pop=="temperate"&tpc$variable=="brate","value"]= tpc[tpc$pop=="temperate"&tpc$variable=="brate","value"]/ max(tpc[tpc$pop=="temperate"&tpc$variable=="brate","value"], na.rm=TRUE)
tpc[tpc$pop=="tropical"&tpc$variable=="brate","value"]= tpc[tpc$pop=="tropical"&tpc$variable=="brate","value"]/ max(tpc[tpc$pop=="tropical"&tpc$variable=="brate","value"], na.rm=TRUE)

#plot
ggplot(data=tpc, aes(x=TempK, y=value, color=variable))+geom_point()+facet_grid(~pop, scales="free")+geom_line()
ylab("fitness component")+xlab("elevation (m)")+geom_smooth(method="loess",se=FALSE, aes(group=year))+theme_bw(base_size = 16)+
  theme(legend.position="bottom")