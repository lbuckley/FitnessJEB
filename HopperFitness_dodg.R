 #Grasshopper biophys
library(ggplot2)
library(maptools)
library(reshape2)
library(ggplot2)
library(MCMCglmm) #for rtnorm function

count=function(x){length(na.omit(x))}

#Model boulderensis at 3048m
#--------------------------
#Read data

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/")
spec.dat=read.csv("SpecData.csv")
#site.list= c("Eldorado","A1","B1","C1")
#specs= c("clav","pell","dodg","sang")
#spec.dat=spec.dat[match(spec.dat$SpecID, specs),]

# Jumping TPC parameters
tpc.dat= read.csv("JumpTPCparams.csv")

#hopping data
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/HopperTPCdata/")
jump.long= read.csv("HoppingData.csv")

specs=c("M. dodgei", "C. pellucida", "M. sanguinipes", "A. clavatus")

#Prefered body temps
Tps= c(32.83, 38.22, 30.63, NA)

#read weather data
#dates are July 5 to September 14
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Grasshoppers/data/PaceDatalogging/Weather2011csv/")
dat=read.csv("DatAll.csv", sep=",", na.strings = c("-49.99","NA"))
#match site data
sites=read.csv("SiteElevations.csv", sep=",")

#combine soil temperatures
dat$SoilTemp= rowMeans(dat[,c("SoilTemp","SoilTemperature")], na.rm = TRUE) 

#add elevation
match1= match(dat$site, sites$Site)
dat$Elev=sites[match1,"Elevation"]

#-----------------------------
#Calculate mass and length, m and g? 
dat$mass_dodg= 1-1.64*10^-4*3048 #*dat$Elev
dat$L_dodg= exp(3.33*0.247*log(dat$mass_dodg))

#source biophys function
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Grasshoppers/biophys/")
source("GrassBiophysFunction_20Jan2012.R")
#biophys=function(Ta, J, Wind, Rad, kt, psi_deg)
#source dtr function 
source("DTRfunction.R")
#source zenith angle calculating function
source("ZenithAngleFunction.R")

#-------------------------
#Fitness functions

#Buckley and Huey SICB
#We estimate fitness as the product of fecundity and survival. 
#Fecundity is quantified as the sum of performance across time steps within a generation, and we assume low but non-zero performance outside the critical thermal limits. 
#We assumed that the probability of survival through a thermal stress event declined exponentially to zero between CTmax and 60 °C. 

#START MORTALITY BEFORE CTmax?
#survival and fecundity  #SURVIVAL COST ABOVE AND BELOW CTmax
#Survival

#Functions

#Need sensitivity analysis
surv<- function(T, CTmin, CTmax, td=4.34){ 
  #10 to 90% of CT range
  CTmin1= CTmin+(CTmax-CTmin)*0.1
  CTmax1= CTmin+(CTmax-CTmin)*0.9
  
  s1= ifelse(T<CTmax1, s<-1, s<- exp(-(T-CTmax1)/td) )
  s2= ifelse(T>CTmin1, s<-1, s<- exp(-(CTmin1-T)/td) )
  s= s1*s2
  return( s*0.8 )
}
plot(0:70, surv(0:70, 10, 60), type="l")

#Model thermoregulation toward Topt
thermoreg.mat<- function(t.mat, Topt){
  Te_sun=t.mat[1]; Te_sh=t.mat[2] 
  To=NA
  if( !is.na(Te_sun) && !is.na(Te_sh)){
    ts= seq(t.mat[2], t.mat[1], 0.1) 
    To= ts[which.min(abs(ts-Topt))]}
  return(To)
}

#Deutsch et al. TPC
#Performance Curve Function from Deutsch et al. 2008
tpc.perf= function(T,Topt,CTmin, CTmax){
  F=T
  F[]=NA
  sigma= (Topt-CTmin)/4
  F[T<=Topt & !is.na(T)]= exp(-((T[T<=Topt & !is.na(T)]-Topt)/(2*sigma))^2) 
  F[T>Topt & !is.na(T)]= 1- ((T[T>Topt & !is.na(T)]-Topt)/(Topt-CTmax))^2
  #1% performance at and outside CT limits
  F[F<=0]<-0.01
  
  return(F)
}


TPC.gausgomp= function(T, To, rho=0.9, sigma, Fmax) Fmax*exp(-exp(rho*(T-To)-6)-sigma*(T-To)^2) # rho and sigma determine, the thermal sensitivity of feeding at temperatures above and below Topt, respectively

#empirical TPC
#CONVERT FROM FT TO M
jump.long$dist =jump.long$dist*0.3048
#aggregate data
dat.pop= aggregate(jump.long, by=list(jump.long$Site, jump.long$Species, jump.long$temp), FUN=mean, na.action = na.omit)  #jump.long$Sex
names(dat.pop)[1:3]= c("Site","Species","temp")
dat.pop1= subset(dat.pop, dat.pop$Species=="dodgei" & dat.pop$elev==3048)

#add CTs to tpc data
dat.tpc= as.data.frame( cbind( c(dat.pop1[,"temp"],7.78,57.16), c(dat.pop1[,"dist"],0,0) ))
colnames(dat.tpc)=c("temp","dist")

lo= loess(dist ~ temp, dat.tpc, span=1)

#---------------
#fitness component plot

par(mfrow=c(1,1))
plot(0:70, surv(0:70, 10, 60), type="l")

#tpc
points(0:70, predict(lo,0:70), type="l", col="green")
points(dat.tpc$temp, dat.tpc$dist)

#================
#Estimate fecundity components

z<-0.001 #specify distance from ground 

#calculate Julian day
date.ct <- strptime(dat$datetime, format="%m/%d/%Y %H:%M", tz="MST")
dat$J= round(as.numeric(julian(date.ct, origin = as.Date("2011-01-01") )))
dat$hour= date.ct$hour

#match site data
match1= match(dat$site, sites$Site) 
dat$lon= sites$Lon[match1]
dat$lat= sites$Lat[match1]

#Calculate zenith
dat$psi=zenith(dat$J, dat$lat, dat$lon, dat$hour)
dat$psi[dat$psi>=85]=85 #set zenith position below horizon to psi=89degrees

#restrict to daylight
dat<- dat[which(dat$hour>5 & dat$hour<17),]

#===============================

#3rd dimension: Te sun, Te shade, Te thermoreg, fecund, surv
#last dimension is individuals
Te.dat= array(data = NA, dim = c(nrow(dat),length(specs),6,500) ) 

spec.k=1

#Estimate Tes, fedund, and surv
    #Estimate Te sun
    #Using soil temperature
    Te.dat[,spec.k,1,1]<-biophys(Ta=dat$SoilTemp, J=dat$J, Wind=dat$Wind, Rad=dat$Rad, kt=1, psi_deg=dat$psi, L=dat$L_dodg, Acondfact=0.0)
    
    #Te IN SHADE
    Te.dat[,spec.k,2,1]<-biophys(Ta=dat$SoilTemp, J=dat$J, Wind=dat$Wind, Rad=0, kt=1, psi_deg=dat$psi, L=dat$L_dodg, Acondfact=0.0)
   
    #Te thermoreg
    Te.dat[,spec.k,3,1]<- apply(cbind(Te.dat[,spec.k,1,1],Te.dat[,spec.k,2,1]), FUN=thermoreg.mat, Topt=spec.dat[spec.k,"PBT"], MARGIN=1)

    #simulate 500 individuals
  inds= which(!is.na(Te.dat[,spec.k,3,1]))
    
  for(ind.k in inds){
  
    #microclimate variability: normal distribution centered at thermoreg temp with sd= (Te_sun-Te_shade)/4, bounded by Te_shade and Te_sun
  ts= rtnorm(n=500, mean = Te.dat[ind.k,spec.k,3,1], sd = 8, lower=(Te.dat[ind.k,spec.k,2,1]-0), upper=(Te.dat[ind.k,spec.k,1,1]+0) )
  #or sd: (Te.dat[ind.k,spec.k,1,1]-Te.dat[ind.k,spec.k,2,1])/4  
  
    #Fecundity
    #Based on CTmin, Topt, CTmax
    Te.dat[ind.k,spec.k,4,]<-tpc.perf( ts, Topt=spec.dat[spec.k,"PBT"], CTmin=spec.dat[spec.k,"CTmin"], CTmax=spec.dat[spec.k,"CTmax"])
    
    #Based on hopping TPC
    # if(!specs[spec.k]=="clav"){ #NOT CLAVATUS
    # ind.elev= dat[ind.k,"Elev"]
    # if(ind.elev==1708)ind.elev<-2195
    # tpc= tpc.dat[which(tpc.dat$spec==specs[spec.k] & tpc.dat$elev_m==2591),]  ##use mid elevation
    # 
    # Te.dat[ind.k,spec.k,6,]<-TPC.gausgomp(ts, To=tpc$To, rho=0.7, sigma=tpc$sigma, Fmax=tpc$Pmax)
    #} #check for clavatus
    Te.dat[ind.k,spec.k,6,]<-predict(lo,ts)
    
    #Survival  
    Te.dat[ind.k,spec.k,5,]<- surv( ts, CTmin= spec.dat[spec.k,"CTmin"], CTmax=spec.dat[spec.k,"CTmax"])  #HEAT STRESS ONLY spec.dat[spec.k,"CTmin"] -20
    
  } #end loops inds
    
#------------------
#ESTIMATE FITNESS
fit= array(NA,dim=c(length(site.list),length(specs),5,500))
#3rd dimension in fecund, surv, fitness

for(site.k in 1:length(site.list) ){
  
  site.inds= which(dat$site==site.list[site.k])
  
  for(spec.k in 1:length(specs) ){

    #estimate fitness as (sum of fecundity)(product of survival)
     fit[site.k, spec.k,1,]= colSums(Te.dat[site.inds,spec.k,4,], na.rm=TRUE) #fecund
     fit[site.k, spec.k,2,]= colSums(Te.dat[site.inds,spec.k,6,], na.rm=TRUE) #fecund tpc
     
      fit[site.k, spec.k,3,]= colMeans(Te.dat[site.inds,spec.k,5,], na.rm=TRUE) #surv

} #end species loop
} #end site loop

#scale fecundity to max
fit[, ,1,]= fit[, ,1,] /max(fit[, ,1,], na.rm=TRUE)
fit[, ,2,]= fit[, ,2,] /max(fit[, ,2,], na.rm=TRUE)

for(site.k in 1:length(site.list) ){
  
  site.inds= which(dat$site==site.list[site.k])
    
    fit[site.k, spec.k,4,]= fit[site.k, spec.k,1,] * fit[site.k, spec.k,3,] #fitness
    fit[site.k, spec.k,5,]= fit[site.k, spec.k,2,] * fit[site.k, spec.k,3,] #fitness tpc
    
} #end site loop

#------------------
#mean fitness across individuals

fit.mean= apply(fit, MARGIN=c(1,2,3), FUN=mean, na.rm=TRUE)

#------------------
#RESHAPE DATA
elevs= c(1708,2195,2591,3048)

f.dat= fit.mean[,1,]
colnames(f.dat)= c("fecundity", "fecundity tpc", "survival", "fitness", "fitness tpc")
f.dat= as.data.frame(f.dat)
f.dat$site= site.list
f.dat$elev= elevs
#melt
f.long= melt(data = f.dat, id.vars = c("site","elev"), measure.vars = c("fecundity","fecundity tpc","survival","fitness", "fitness tpc"))
#order fitness factors
f.long$component<- f.long$variable
f.long$component= factor(f.long$component, levels=c("fecundity","fecundity tpc","survival","fitness", "fitness tpc") )

#PLOT
#selections
f.long= f.long[f.long$component %in% c("fecundity tpc","survival","fitness tpc"),]

fit.plot= ggplot(data=f.long, aes(x=elev, y = value, color=component))+geom_line()+
  theme_bw()+theme(legend.position="bottom") +ylab("fitness component")+xlab("elevation (m)")
#+scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
#                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3"))

#FIGURE
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/figures/")
pdf("FitnessFig.pdf", height = 10, width = 6)
fit.plot
dev.off()
