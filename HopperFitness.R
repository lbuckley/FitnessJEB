 #Grasshopper biophys
library(ggplot2)
library(maptools)

count=function(x){length(na.omit(x))}

#--------------------------
#Grasshopper data

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/")
spec.dat=read.csv("SpecData.csv")

#------------------------------

#read weather data
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

# Function to plot line if trend is significant
PlotSigTrend=function(x,y, col1, lty1){
#mod is linear model

if( sum(!is.na(y))>0){ #check data exists
mod<-lm(y~x)
mods<-summary(mod)
f <- mods$fstatistic
p<-pf(f[1], f[2], f[3], lower=FALSE) 
if(!is.na(p))if(p<0.05)abline(mod, col=col1, lty=lty1)
}}

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

#------------------------
#Make 1/kT column
k=8.62*10^-5 #eVK^-1
dat$kT= 1/(k* (dat$Temp+273))

#Match lengths and mass across sites
dat$L_clav=rep(0, nrow(dat))
dat$L_pell=rep(0, nrow(dat))
dat$L_dodg=rep(0, nrow(dat))
dat$L_sang=rep(0, nrow(dat))
dat$mass_clav=rep(0, nrow(dat))
dat$mass_pell=rep(0, nrow(dat))
dat$mass_dodg=rep(0, nrow(dat))
dat$mass_sang=rep(0, nrow(dat))

#Eldorado
which.site= which(dat$site=="Eldorado")
dat$L_clav[which.site]=spec.dat[which(spec.dat$SpecID=="clav"),"Lmm_Eldo"]/1000
dat$L_pell[which.site]=spec.dat[which(spec.dat$SpecID=="pell"),"Lmm_Eldo"]/1000
dat$L_dodg[which.site]=spec.dat[which(spec.dat$SpecID=="dodg"),"Lmm_Eldo"]/1000
dat$L_sang[which.site]=spec.dat[which(spec.dat$SpecID=="sang"),"Lmm_Eldo"]/1000
dat$mass_clav[which.site]=spec.dat[which(spec.dat$SpecID=="clav"),"Massg_Eldo"]/1000
dat$mass_pell[which.site]=spec.dat[which(spec.dat$SpecID=="pell"),"Massg_Eldo"]/1000
dat$mass_dodg[which.site]=spec.dat[which(spec.dat$SpecID=="dodg"),"Massg_Eldo"]/1000
dat$mass_sang[which.site]=spec.dat[which(spec.dat$SpecID=="sang"),"Massg_Eldo"]/1000

#A1
which.site= which(dat$site=="A1")
dat$L_clav[which.site]=spec.dat[which(spec.dat$SpecID=="clav"),"Lmm_A1"]/1000
dat$L_pell[which.site]=spec.dat[which(spec.dat$SpecID=="pell"),"Lmm_A1"]/1000
dat$L_dodg[which.site]=spec.dat[which(spec.dat$SpecID=="dodg"),"Lmm_A1"]/1000
dat$L_sang[which.site]=spec.dat[which(spec.dat$SpecID=="sang"),"Lmm_A1"]/1000
dat$mass_clav[which.site]=spec.dat[which(spec.dat$SpecID=="clav"),"Massg_A1"]/1000
dat$mass_pell[which.site]=spec.dat[which(spec.dat$SpecID=="pell"),"Massg_A1"]/1000
dat$mass_dodg[which.site]=spec.dat[which(spec.dat$SpecID=="dodg"),"Massg_A1"]/1000
dat$mass_sang[which.site]=spec.dat[which(spec.dat$SpecID=="sang"),"Massg_A1"]/1000

#B1
which.site= which(dat$site=="B1")
dat$L_clav[which.site]=spec.dat[which(spec.dat$SpecID=="clav"),"Lmm_B1"]/1000
dat$L_pell[which.site]=spec.dat[which(spec.dat$SpecID=="pell"),"Lmm_B1"]/1000
dat$L_dodg[which.site]=spec.dat[which(spec.dat$SpecID=="dodg"),"Lmm_B1"]/1000
dat$L_sang[which.site]=spec.dat[which(spec.dat$SpecID=="sang"),"Lmm_B1"]/1000
dat$mass_clav[which.site]=spec.dat[which(spec.dat$SpecID=="clav"),"Massg_B1"]/1000
dat$mass_pell[which.site]=spec.dat[which(spec.dat$SpecID=="pell"),"Massg_B1"]/1000
dat$mass_dodg[which.site]=spec.dat[which(spec.dat$SpecID=="dodg"),"Massg_B1"]/1000
dat$mass_sang[which.site]=spec.dat[which(spec.dat$SpecID=="sang"),"Massg_B1"]/1000

#C1
which.site= which(dat$site=="C1")
dat$L_clav[which.site]=spec.dat[which(spec.dat$SpecID=="clav"),"Lmm_C1"]/1000
dat$L_pell[which.site]=spec.dat[which(spec.dat$SpecID=="pell"),"Lmm_C1"]/1000
dat$L_dodg[which.site]=spec.dat[which(spec.dat$SpecID=="dodg"),"Lmm_C1"]/1000
dat$L_sang[which.site]=spec.dat[which(spec.dat$SpecID=="sang"),"Lmm_C1"]/1000
dat$mass_clav[which.site]=spec.dat[which(spec.dat$SpecID=="clav"),"Massg_C1"]/1000
dat$mass_pell[which.site]=spec.dat[which(spec.dat$SpecID=="pell"),"Massg_C1"]/1000
dat$mass_dodg[which.site]=spec.dat[which(spec.dat$SpecID=="dodg"),"Massg_C1"]/1000
dat$mass_sang[which.site]=spec.dat[which(spec.dat$SpecID=="sang"),"Massg_C1"]/1000

#-----------------------------
#Calculate mass and length 
dat$mass_clav= 0.42-8.15*10^-5*dat$Elev
dat$mass_pell= 0.77-1.48*10^-4*dat$Elev
dat$mass_dodg= 1-1.64*10^-4*dat$Elev
dat$mass_sang= 0.456-2.74*10^-5*dat$Elev
dat$L_clav=  exp(3.33*0.247*log(dat$mass_clav)) #in mm
dat$L_pell= exp(3.33*0.247*log(dat$mass_pell))
dat$L_dodg= exp(3.33*0.247*log(dat$mass_dodg))
dat$L_sang= exp(3.33*0.247*log(dat$mass_sang))

####################################

##Calculate Te for each species
#Using soil temperature
dat$Te_clav<-biophys(Ta=dat$SoilTemp, J=dat$J, Wind=dat$Wind, Rad=dat$Rad, kt=1, psi_deg=dat$psi, L=dat$L_clav, Acondfact=0.0)
dat$Te_pell<-biophys(dat$SoilTemp, dat$J, dat$Wind, dat$Rad, kt=1, dat$psi, L=dat$L_pell, Acondfact=0.0)
dat$Te_dodg<-biophys(dat$SoilTemp, dat$J, dat$Wind, dat$Rad, kt=1, dat$psi, L=dat$L_dodg, Acondfact=0.0)
dat$Te_sang<-biophys(dat$SoilTemp, dat$J, dat$Wind, dat$Rad, kt=1, dat$psi, L=dat$L_sang, Acondfact=0.0)

dat$Te_clav[is.na(dat$Te_clav)]=NA
dat$Te_pell[is.na(dat$Te_pell)]=NA
dat$Te_dodg[is.na(dat$Te_dodg)]=NA
dat$Te_sang[is.na(dat$Te_sang)]=NA

#Te IN SHADE
dat$Te_clav.shade<-biophys(Ta=dat$SoilTemp, J=dat$J, Wind=dat$Wind, Rad=0, kt=1, psi_deg=dat$psi, L=dat$L_clav, Acondfact=0.0)
dat$Te_pell.shade<-biophys(dat$SoilTemp, dat$J, dat$Wind, Rad=0, kt=1, dat$psi, L=dat$L_pell, Acondfact=0.0)
dat$Te_dodg.shade<-biophys(dat$SoilTemp, dat$J, dat$Wind, Rad=0, kt=1, dat$psi, L=dat$L_dodg, Acondfact=0.0)
dat$Te_sang.shade<-biophys(dat$SoilTemp, dat$J, dat$Wind, Rad=0, kt=1, dat$psi, L=dat$L_sang, Acondfact=0.0)

dat$Te_clav.shade[is.na(dat$Te_clav.shade)]=NA
dat$Te_pell.shade[is.na(dat$Te_pell.shade)]=NA
dat$Te_dodg.shade[is.na(dat$Te_dodg.shade)]=NA
dat$Te_sang.shade[is.na(dat$Te_sang.shade)]=NA

#=====================================
#Biological calculations

##Calculate metabolic rate, VCO2mLh
#clav
#wrow=which(spec.dat$SpecID=="clav")
#dat$MR_clav= exp( spec.dat[wrow,"b0"]+spec.dat[wrow,"b1"]*log(dat$mass_clav)+spec.dat[wrow,"b2"]*(1/(k* (dat$Te_clav+273)))+spec.dat[wrow,"b3"]*dat$Elev)

##Calculate activity time
#dat$AT_clav[which(dat$Te_clav >= spec.dat[wrow, "CTmin"] & dat$Te_clav <= spec.dat[wrow, "CTmax"])]=1
#dat$AT_clav[which(dat$Te_clav >= spec.dat[wrow, "Tb20"] & dat$Te_clav <= spec.dat[wrow, "Tb80"])]=1

#Estimate Te

#Estiamte survival

#Estimate fitness




