 #Grasshopper biophys
library(ggplot2)
library(maptools)
library(reshape2)
library(ggplot2)
library(MCMCglmm) #for rtnorm function

count=function(x){length(na.omit(x))}

#-------------------------
#Fitness functions

#Buckley and Huey SICB
#We estimate fitness as the product of fecundity and survival. 
#Fecundity is quantified as the sum of performance across time steps within a generation, and we assume low but non-zero performance outside the critical thermal limits. 
#We assumed that the probability of survival through a thermal stress event declined exponentially to zero between CTmax and 60 °C. 

#START MORTALITY BEFORE CTmax?
#survival and fecundity  #SURVIVAL COST ABOVE AND BELOW CTmax
#Survival

#Need sensitivity analysis
surv<- function(T, CTmin, CTmax, td=4.34){ 
  #10 to 90% of CT range
  CTmin1= CTmin+(CTmax-CTmin)*0.2
  CTmax1= CTmin+(CTmax-CTmin)*0.8
  
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

#--------------------------
#Grasshopper data

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/")
spec.dat=read.csv("SpecData.csv")

# Jumping TPC parameters
tpc.dat= read.csv("JumpTPCparams.csv")
# or JumpModels_2019

#PLOT DATA
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/figures/")
pdf("TPCFig.pdf", height = 6, width = 10)

specs1= c("boulderensis","pellucida","sanguinipes")
elevs1=c(2195,2591,3048)
temps=0:60

elev_lab= paste(elevs, "m", sep="")
specs_lab=c("M. boulderensis", "C. pellucida", "M. sanguinipes", "A. clavatus")
par(mfrow=c(1,3), cex=1.2, lwd=2, mar=c(1, 2.5, 2.5, 0.2), mgp=c(1.5, 0.5, 0), oma=c(2,2,2,0), bty="l")

for(sp.k in 1:3){
for(elev.k in 1:3){
  
  tpc= tpc.dat[which(tpc.dat$species==specs1[sp.k] & tpc.dat$elev_m==elevs1[elev.k]),]
  if(elev.k==1) plot(temps, TPC.gausgomp(temps, To=tpc$To, rho=0.7, sigma=tpc$sigma, Fmax=tpc$Pmax), type="l", col="red", ylim=range(0,0.6),xlim=range(0,60), ylab="Hopping distance (cm)", main=specs1[sp.k] )
  if(elev.k==2) points(temps, TPC.gausgomp(temps, To=tpc$To, rho=0.7, sigma=tpc$sigma, Fmax=tpc$Pmax), type="l", col="green")
  if(elev.k==3) points(temps, TPC.gausgomp(temps, To=tpc$To, rho=0.7, sigma=tpc$sigma, Fmax=tpc$Pmax), type="l", col="blue")
  
} #end sites
  
#add CTmin and CTmax  
points(c(spec.dat[sp.k,"CTmin"],spec.dat[sp.k,"PBT"],spec.dat[sp.k,"CTmax"]  ), c(0,0,0))
    
} #end species

#add legend
legend(35, 0.6, legend=c("2195m", "2591m", "3048m"), col=c("red", "green", "blue"), lty=1, cex=0.8)

dev.off()

#------------------------------

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

#source biophys function
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Grasshoppers/biophys/")
source("GrassBiophysFunction_20Jan2012.R")
#biophys=function(Ta, J, Wind, Rad, kt, psi_deg)
#source dtr function 
source("DTRfunction.R")
#source zenith angle calculating function
source("ZenithAngleFunction.R")

#------------------------
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

#-----------------------------
#Calculate mass and length, m and g? 
dat$mass_clav= 0.42-8.15*10^-5*2591 #*dat$Elev
dat$mass_pell= 0.77-1.48*10^-4*2591 #*dat$Elev
dat$mass_dodg= 1-1.64*10^-4*2591 #*dat$Elev
dat$mass_sang= 0.456-2.74*10^-5*2591 #*dat$Elev
dat$L_clav=  exp(3.33*0.247*log(dat$mass_clav)) #in mm
dat$L_pell= exp(3.33*0.247*log(dat$mass_pell))
dat$L_dodg= exp(3.33*0.247*log(dat$mass_dodg))
dat$L_sang= exp(3.33*0.247*log(dat$mass_sang))

#===============================

site.list= c("Eldorado","A1","B1","C1")
specs= c("clav","pell","dodg","sang")
spec.dat=spec.dat[match(spec.dat$SpecID, specs),]

#3rd dimension: Te sun, Te shade, Te thermoreg, fecund, surv
#last dimension is individuals
Te.dat= array(data = NA, dim = c(nrow(dat),length(specs),6,500) ) 

#Estimate Tes, fedund, and surv
    for(spec.k in 1:length(specs) ){
  
    #Estimate Te sun
    #Using soil temperature
    Te.dat[,spec.k,1,1]<-biophys(Ta=dat$SoilTemp, J=dat$J, Wind=dat$Wind, Rad=dat$Rad, kt=1, psi_deg=dat$psi, L=dat[,(25+spec.k)], Acondfact=0.0)
    
    #Te IN SHADE
    Te.dat[,spec.k,2,1]<-biophys(Ta=dat$SoilTemp, J=dat$J, Wind=dat$Wind, Rad=0, kt=1, psi_deg=dat$psi, L=dat[,(25+spec.k)], Acondfact=0.0)
   
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
    if(!specs[spec.k]=="clav"){ #NOT CLAVATUS
    ind.elev= dat[ind.k,"Elev"]
    if(ind.elev==1708)ind.elev<-2195
    tpc= tpc.dat[which(tpc.dat$spec==specs[spec.k] & tpc.dat$elev_m==2591),]  ##use mid elevation
    
    Te.dat[ind.k,spec.k,6,]<-TPC.gausgomp(ts, To=tpc$To, rho=0.7, sigma=tpc$sigma, Fmax=tpc$Pmax)
    } #check for clavatus
    
    #Survival  
    Te.dat[ind.k,spec.k,5,]<- surv( ts, CTmin= spec.dat[spec.k,"CTmin"], CTmax=spec.dat[spec.k,"CTmax"])  #HEAT STRESS ONLY spec.dat[spec.k,"CTmin"] -20
    
  } #end loops time points
    
    } #end spec loop  

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
  
  for(spec.k in 1:length(specs) ){
    
    fit[site.k, spec.k,4,]= fit[site.k, spec.k,1,] * fit[site.k, spec.k,3,] #fitness
    fit[site.k, spec.k,5,]= fit[site.k, spec.k,2,] * fit[site.k, spec.k,3,] #fitness tpc
    
  } #end species loop
} #end site loop

#------------------
#mean fitness across individuals

fit.mean= apply(fit, MARGIN=c(1,2,3), FUN=mean, na.rm=TRUE)

#------------------
#RESHAPE DATA
elevs= c(1708,2195,2591,3048)

f.dat= rbind(fit.mean[,,1],fit.mean[,,2],fit.mean[,,3],fit.mean[,,4],fit.mean[,,5])
colnames(f.dat)= specs
f.dat= as.data.frame(f.dat)
f.dat$component= c(rep("fecundity", length(site.list)),rep("fecundity tpc", length(site.list)), rep("survival", length(site.list)), rep("fitness", length(site.list)), rep("fitness tpc", length(site.list)) )
f.dat$site= rep(site.list, 5)
f.dat$elev= rep(elevs, 5)
#melt
f.long= melt(data = f.dat, id.vars = c("component","site","elev"), measure.vars = c("clav", "pell", "dodg", "sang"))
colnames(f.long)[4]<-"species"
#order fitness factors
f.long$component= factor(f.long$component, levels=c("fecundity","fecundity tpc","survival","fitness", "fitness tpc") )

#PLOT
#selections
f.long= f.long[f.long$species %in% c("pell", "dodg", "sang"),]
f.long= f.long[f.long$component %in% c("fecundity tpc","survival","fitness tpc"),]

fit.plot= ggplot(data=f.long, aes(x=elev, y = value, color=component))+geom_line()+facet_wrap(~species, ncol=1, scales="free_y")+
  theme_bw()+theme(legend.position="bottom") +ylab("fitness component")+xlab("elevation (m)")
#+scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
#                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3"))

#FIGURE
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/figures/")
pdf("FitnessFig.pdf", height = 10, width = 6)
fit.plot
dev.off()

#=====================================
#Biological calculations

##Calculate metabolic rate, VCO2mLh
#clav
#wrow=which(spec.dat$SpecID=="clav")
#dat$MR_clav= exp( spec.dat[wrow,"b0"]+spec.dat[wrow,"b1"]*log(dat$mass_clav)+spec.dat[wrow,"b2"]*(1/(k* (dat$Te_clav+273)))+spec.dat[wrow,"b3"]*dat$Elev)

##Calculate activity time
#dat$AT_clav[which(dat$Te_clav >= spec.dat[wrow, "CTmin"] & dat$Te_clav <= spec.dat[wrow, "CTmax"])]=1
#dat$AT_clav[which(dat$Te_clav >= spec.dat[wrow, "Tb20"] & dat$Te_clav <= spec.dat[wrow, "Tb80"])]=1

#=====================================
#BOULDERENSIS
#Make new boulderensis TPC, 3048m?

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/HopperTPCdata/")
hop= read.csv("HoppingData.csv")
jump.long=hop



specs1= c("boulderensis","pellucida","sanguinipes")
elevs1=c(2195,2591,3048)
sp.k=2
elev.k=2
temps=0:60

#estimate TPC
tpc= tpc.dat[which(tpc.dat$species==specs1[sp.k] & tpc.dat$elev_m==elevs1[elev.k]),]
perf= 

par(mfrow=c(1,1), cex=1.2, lwd=2, mar=c(1, 2.5, 2.5, 0.2), mgp=c(1.5, 0.5, 0), oma=c(2,2,2,0), bty="l")
TPC.plot= 
  

    
    if(elev.k==1) plot(temps, TPC.gausgomp(temps, To=tpc$To, rho=0.7, sigma=tpc$sigma, Fmax=tpc$Pmax), type="l", col="red", ylim=range(0,0.6),xlim=range(0,60), ylab="Hopping distance (cm)", main=specs1[sp.k] )
    if(elev.k==2) points(temps, TPC.gausgomp(temps, To=tpc$To, rho=0.7, sigma=tpc$sigma, Fmax=tpc$Pmax), type="l", col="green")
    if(elev.k==3) points(temps, TPC.gausgomp(temps, To=tpc$To, rho=0.7, sigma=tpc$sigma, Fmax=tpc$Pmax), type="l", col="blue")
    
  } #end sites
  
  #add CTmin and CTmax  
  points(c(spec.dat[sp.k,"CTmin"],spec.dat[sp.k,"PBT"],spec.dat[sp.k,"CTmax"]  ), c(0,0,0))
  
} #end species


#*********************************
#ANALYSIS OF THE TEMPERATURE DEPENDENCE OF PERFORMANCE AND FEEDING RATES

#$$$$$$ THERMOREG TOGLE
TR=TRUE
#$$$$$$$$

library(car)
library(nls2) #for TPC fitting
library(minpack.lm)

count=function(x) length(x)

#HOPPING
sites= c("Redfox", "A1", "B1", "C1", "D1")  
elevs= c(1574, 2195, 2591, 3048, 3739)

#read in data
setwd("F:\\Work\\HopperTPC\\Data\\DataClean\\")
jump.long= read.csv("HoppingData.csv")

specs=c("M. dodgei", "C. pellucida", "M. sanguinipes", "A. clavatus")

#Prefered body temps
Tps= c(32.83, 38.22, 30.63, NA)

#--------------------------
#CONVERT FROM FT TO M
jump.long$dist =jump.long$dist*0.3048

#aggregate long data
dat.pop= aggregate(jump.long, by=list(jump.long$Site, jump.long$Species, jump.long$temp), FUN=mean, na.action = na.omit)  #jump.long$Sex
dat.pop.sd= aggregate(jump.long, by=list(jump.long$Site, jump.long$Species, jump.long$temp), FUN=sd, na.rm = TRUE)
dat.pop.count= aggregate(jump.long, by=list(jump.long$Site, jump.long$Species, jump.long$temp), FUN=count)
names(dat.pop)[1:3]= c("Site","Species","temp")
dat.pop.se= dat.pop.count
dat.pop.se[,9:12]= dat.pop.sd[,9:12]/sqrt(dat.pop.count[,9:12])
names(dat.pop.se)[1:3]= c("Site","Species","temp")

#FIGURE 2
elev_lab= paste(elevs, "m", sep="")
specs= c("dodgei", "pellucida", "sanguinipes", "clavatus")
specs_lab=c("M. boulderensis \ncool adapted", "C. pellucida \nwarm adapted", "M. sanguinipes \ngeneralist", "A. clavatus")
par(mfrow=c(1,3), cex=1.4, lwd=2, mar=c(1, 1, 1.6, 0.7), mgp=c(1.5, 0.5, 0), oma=c(2,2,0,0), bty="l")

cols=c("darkgreen", "limegreen", "turquoise", "dodgerblue", "darkblue")
#colfunc <- colorRampPalette(c("green", "blue"))
#cols=colfunc(5)
pch.match=c("f","m")
pchs=c(1,19) 
#pch=pchs[match(dat.spec1$Sex, pch.match)]

for(i in 1:3){
  dat.spec=subset(dat.pop, dat.pop$Species==specs[i])
  dat.spec.se=subset(dat.pop.se, dat.pop.se$Species==specs[i])
  
  for(site in 1:5){
    dat.spec1=subset(dat.spec, dat.spec$Site==sites[site])
    dat.spec1.se=subset(dat.spec.se, dat.spec.se$Site==sites[site])
    
    if(site==1){
      plot(dat.spec1$temp, dat.spec1$dist, col=cols[site], xlim=range(10,35), ylim=range(0.13,0.63), xlab="", ylab="", pch=1, main="", type="b") 
      title(specs_lab[i], line=-0.9)
    }
    if(site>1) points(dat.spec1$temp, dat.spec1$dist, col=cols[site], type="b", pch=1)
    
    arrows(dat.spec1$temp, dat.spec1$dist-dat.spec1.se$dist, dat.spec1$temp, dat.spec1$dist+dat.spec1.se$dist, code=3, angle=90,length=0.1, col=cols[site])
  } #end loop sites
  
  if(i==1) legend("bottomright", elev_lab, pch=1, col=cols,ncol=2, bty="n")
}

mtext("Distance (m)", side = 2, line = 0.5, outer = TRUE, cex=1.5)
mtext("Temperature (?C)", side = 1, line = 0.5, outer = TRUE, cex=1.5)


