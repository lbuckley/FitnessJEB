 #Grasshopper biophys
library(ggplot2)
library(maptools)
library(reshape2)
library(ggplot2)
library(MCMCglmm) #for rtnorm function

count=function(x){length(na.omit(x))}

#--------------------------
#Grasshopper data

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/")
spec.dat=read.csv("SpecData.csv")

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

#-----------------------------
#Calculate mass and length, m and g? 
dat$mass_clav= 0.42-8.15*10^-5*dat$Elev
dat$mass_pell= 0.77-1.48*10^-4*dat$Elev
dat$mass_dodg= 1-1.64*10^-4*dat$Elev
dat$mass_sang= 0.456-2.74*10^-5*dat$Elev
dat$L_clav=  exp(3.33*0.247*log(dat$mass_clav)) #in mm
dat$L_pell= exp(3.33*0.247*log(dat$mass_pell))
dat$L_dodg= exp(3.33*0.247*log(dat$mass_dodg))
dat$L_sang= exp(3.33*0.247*log(dat$mass_sang))

#===============================
#Fitness functions

#Buckley and Huey SICB
#We estimate fitness as the product of fecundity and survival. 
#Fecundity is quantified as the sum of performance across time steps within a generation, and we assume low but non-zero performance outside the critical thermal limits. 
#We assumed that the probability of survival through a thermal stress event declined exponentially to zero between CTmax and 60 °C. 


#START MORTALITY BEFORE CTmax?
#survival and fecundity  #SURVIVAL COST ABOVE AND BELOW CTmax
#Survival

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

#-------------------
site.list= c("Eldorado","A1","B1","C1")
specs= c("clav","pell","dodg","sang")
spec.dat=spec.dat[match(spec.dat$SpecID, specs),]

#3rd dimension: Te sun, Te shade, Te thermoreg, fecund, surv
#last dimension is individuals
Te.dat= array(data = NA, dim = c(nrow(dat),length(specs),5,500) ) 

#Estimate Tes, fedund, and surv
    for(spec.k in 1:length(specs) ){
  
    #Estimate Te sun
    #Using soil temperature
    Te.dat[,spec.k,1,1]<-biophys(Ta=dat$SoilTemp, J=dat$J, Wind=dat$Wind, Rad=dat$Rad, kt=1, psi_deg=dat$psi, L=dat[,(25+spec.k)], Acondfact=0.0)
    
    #Te IN SHADE
    Te.dat[,spec.k,2,1]<-biophys(Ta=dat$SoilTemp, J=dat$J, Wind=dat$Wind, Rad=0, kt=1, psi_deg=dat$psi, L=dat[,(25+spec.k)], Acondfact=0.0)
   
    #Te thermoreg
    Te.dat[,spec.k,3,1]<-apply(cbind(Te.dat[,spec.k,1,1],Te.dat[,spec.k,2,1]), FUN=thermoreg.mat, Topt=spec.dat[spec.k,"PBT"], MARGIN=1)

    #simulate 500 individuals
  inds= which(!is.na(Te.dat[,spec.k,3,1]))
    
  for(ind.k in inds){
  
    #microcliamte variability: normal distribution centered at thermoreg temp with sd= (Te_sun-Te_shade)/4, bounded by Te_shade -5 and Te_sun +5
  ts= rtnorm(n=500, mean = Te.dat[ind.k,spec.k,3,1], sd = 4, lower=(Te.dat[ind.k,spec.k,2,1]-0), upper=(Te.dat[ind.k,spec.k,1,1]+0) )
  #or sd: (Te.dat[ind.k,spec.k,1,1]-Te.dat[ind.k,spec.k,2,1])/4  
  
    #Fecundity
    Te.dat[ind.k,spec.k,4,]<-tpc.perf( ts, Topt=spec.dat[spec.k,"PBT"], CTmin=spec.dat[spec.k,"CTmin"], CTmax=spec.dat[spec.k,"CTmax"])
    
    #Survival  
    Te.dat[ind.k,spec.k,5,]<- surv( ts, CTmin= -20, CTmax=spec.dat[spec.k,"CTmax"])  #HEAT STRESS ONLY spec.dat[spec.k,"CTmin"]
    
  } #end loops time points
    
    } #end spec loop  

#------------------
#ESTIMATE FITNESS
fit= array(NA,dim=c(length(site.list),length(specs),3,500))
#3rd dimension in fecund, surv, fitness

for(site.k in 1:length(site.list) ){
  
  site.inds= which(dat$site==site.list[site.k])
  
  for(spec.k in 1:length(specs) ){

    #estimate fitness as (sum of fecundity)(product of survival)
     fit[site.k, spec.k,1,]= colSums(Te.dat[site.inds,spec.k,4,], na.rm=TRUE) #fecund
      fit[site.k, spec.k,2,]= colMeans(Te.dat[site.inds,spec.k,5,], na.rm=TRUE) #surv
      fit[site.k, spec.k,3,]= fit[site.k, spec.k,1,] * fit[site.k, spec.k,2,] #fitness
} #end species loop
} #end site loop

#------------------
#mean fitness across individuals

fit.mean= apply(fit, MARGIN=c(1,2,3), FUN=mean, na.rm=TRUE)

#------------------
#RESHAPE DATA
elevs= c(1708,2195,2591,3048)

f.dat= rbind(fit.mean[,,1],fit.mean[,,2],fit.mean[,,3])
colnames(f.dat)= specs
f.dat= as.data.frame(f.dat)
f.dat$component= c(rep("fecundity", length(site.list)), rep("survival", length(site.list)), rep("fitness", length(site.list)) )
f.dat$site= rep(site.list, 3)
f.dat$elev= rep(elevs, 3)
#melt
f.long= melt(data = f.dat, id.vars = c("component","site","elev"), measure.vars = c("clav", "pell", "dodg", "sang"))
colnames(f.long)[4]<-"species"
#order fitness factors
f.long$component= factor(f.long$component, levels=c("fecundity","survival","fitness") )

#PLOT
fit.plot= ggplot(data=f.long, aes(x=elev, y = value, color=species))+geom_line()+facet_wrap(~component, ncol=1, scales="free_y")+
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






