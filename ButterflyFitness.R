library(ggplot2)

#CHOOSE PROJECTION
proj.k=2 #1: bcc-csm1-1.1.rcp60, 2: ccsm4.1.rcp60, 3: gfdl-cm3.1.rcp60
projs=c("bcc-csm","ccsm4","gfdl")

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/ColiasBiogeog/OUT/")

#read data
filename= paste("lambda1_",projs[proj.k],".rds",sep="")
lambda <- readRDS(filename)
#Lambda[yr.k, cell.k, abs.k, gen.k, ]= c(Lambda1, mean(FAT.ind), mean(ev.ind), mean(T.ind, na.rm=T) )
#ev.ind: calculate geometric mean of egg viability within flight period
#FAT.ind: average FAT over days

Tp= readRDS(paste("PupTemps_",projs[proj.k],".rds",sep=""))
pts= read.csv(paste("COpoints_",projs[proj.k],".rds",sep=""))

#--------------------
#For given time period, plot FAT and ev.ind across elevation

fec= lambda[150, , 3, 2, 2]
surv= lambda[150, , 3, 2, 3]
dat= cbind(pts, fec, surv)

ggplot(data=dat, aes(x=elev, y=fec))+geom_point()
ggplot(data=dat, aes(x=elev, y=surv))+geom_point()

#PLOT OUT TWO TIME PERIODS
#TWO ABSORPTIVITIES?
#GENERATIONS?

