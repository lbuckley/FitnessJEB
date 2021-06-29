library(ggplot2)
library(cowplot)
library(reshape2)

#Prepare Colias lambda data

years= 1950:2099
#absorptivity
aseq= seq(0.4,0.7,0.05)

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

fec= lambda[50, , 3, 2, 2]
surv= lambda[50, , 3, 2, 3]
dat1= cbind(pts, fec, surv)
dat1$year= 50
dat1$abs=3

fec= lambda[150, , 3, 2, 2]
surv= lambda[150, , 3, 2, 3]
dat2= cbind(pts, fec, surv)
dat2$year= 150
dat2$abs=3

fec= lambda[50, , 6, 2, 2]
surv= lambda[50, , 6, 2, 3]
dat3= cbind(pts, fec, surv)
dat3$year= 50
dat3$abs=6

fec= lambda[150, , 6, 2, 2]
surv= lambda[150, , 6, 2, 3]
dat4= cbind(pts, fec, surv)
dat4$year= 150
dat4$abs=6

dat=rbind(dat1, dat2, dat3, dat4)
#combine surv and fecundity
ldat= melt(dat, id.vars = c("lon","lat","lon.ind","lat.ind","ind","elev","airpr","year","abs"), measure.vars = c("fec", "surv"))
ldat=na.omit(ldat)
ldat$absorp= paste("wing absorptivity=", aseq[ldat$abs])
ldat$fitcomp= "flight activity time (hr)"
ldat$fitcomp[ldat$variable=="surv"]= "egg viability (%)"
ldat$fitcomp= ordered(ldat$fitcomp, levels=c("flight activity time (hr)","egg viability (%)") )
ldat$year= years[ldat$year]
ldat$year= as.factor(ldat$year)

write.csv(ldat,"ColiasLambdaData.csv")