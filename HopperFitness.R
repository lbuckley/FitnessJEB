 #Grasshopper biophys
library(maptools)

count=function(x){length(na.omit(x))}

#--------------------------
#Grasshopper data

setwd("\\\\Bioark.bio.unc.edu\\buckleylauren\\Work\\Grasshoppers\\data\\MRdata\\Data_fin\\")
spec.dat=read.csv("SpecData.csv")

#------------------------------

#read weather data
setwd("\\\\Bioark.bio.unc.edu\\buckleylauren\\Work\\Grasshoppers\\data\\PaceDatalogging\\Weather2011csv\\")
dat=read.csv("DatAll.csv", sep=",", na.strings = c("-49.99","NA"))
#match site data
sites=read.csv("SiteElevations.csv", sep=",")
 
#combine soil temperatures
dat$SoilTemp= rowMeans(dat[,c("SoilTemp","SoilTemperature")], na.rm = TRUE) 

#add elevation
match1= match(dat$site, sites$Site)
dat$Elev=sites[match1,"Elevation"]

#source biophys function
setwd("\\\\Bioark\\buckleylauren\\Work\\Grasshoppers\\biophys\\")
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

#Test plot
#plot(dat$SoilTemp, dat$TeSoil, ylim=range(0,60), xlim=range(0,60))
#abline(coef=c(0,1), col="red")

#-----------------------------
##Calculate metabolic rate, VCO2mLh

#clav
wrow=which(spec.dat$SpecID=="clav")
dat$MR_clav= exp( spec.dat[wrow,"b0"]+spec.dat[wrow,"b1"]*log(dat$mass_clav)+spec.dat[wrow,"b2"]*(1/(k* (dat$Te_clav+273)))+spec.dat[wrow,"b3"]*dat$Elev)

#pell
wrow=which(spec.dat$SpecID=="pell")
dat$MR_pell= exp( spec.dat[wrow,"b0"]+spec.dat[wrow,"b1"]*log(dat$mass_pell)+spec.dat[wrow,"b2"]*(1/(k* (dat$Te_pell+273)))+spec.dat[wrow,"b3"]*dat$Elev)

#dodg
wrow=which(spec.dat$SpecID=="dodg")
dat$MR_dodg= exp( spec.dat[wrow,"b0"]+spec.dat[wrow,"b1"]*log(dat$mass_dodg)+spec.dat[wrow,"b2"]*(1/(k* (dat$Te_dodg+273)))+spec.dat[wrow,"b3"]*dat$Elev)

#sang
wrow=which(spec.dat$SpecID=="sang")
dat$MR_sang= exp( spec.dat[wrow,"b0"]+spec.dat[wrow,"b1"]*log(dat$mass_sang)+spec.dat[wrow,"b2"]*(1/(k* (dat$Te_sang+273)))+spec.dat[wrow,"b3"]*dat$Elev)

#-----------------------------
##Calculate activity time
dat$AT_clav=rep(0, nrow(dat))
dat$AT_pell=rep(0, nrow(dat))
dat$AT_dodg=rep(0, nrow(dat))
dat$AT_sang=rep(0, nrow(dat))

#clav
wrow=which(spec.dat$SpecID=="clav")
dat$AT_clav[which(dat$Te_clav >= spec.dat[wrow, "CTmin"] & dat$Te_clav <= spec.dat[wrow, "CTmax"])]=1
#pell
wrow=which(spec.dat$SpecID=="pell")
dat$AT_pell[which(dat$Te_pell >= spec.dat[wrow, "CTmin"] & dat$Te_pell <= spec.dat[wrow, "CTmax"])]=1
#dodg
wrow=which(spec.dat$SpecID=="dodg")
dat$AT_dodg[which(dat$Te_dodg >= spec.dat[wrow, "CTmin"] & dat$Te_dodg <= spec.dat[wrow, "CTmax"])]=1
#sang
wrow=which(spec.dat$SpecID=="sang")
dat$AT_sang[which(dat$Te_sang >= spec.dat[wrow, "CTmin"] & dat$Te_sang <= spec.dat[wrow, "CTmax"])]=1

#USE 20% and 80% quantiles
dat$AT_clav=rep(0, nrow(dat))
dat$AT_pell=rep(0, nrow(dat))
dat$AT_dodg=rep(0, nrow(dat))
dat$AT_sang=rep(0, nrow(dat))

#clav
wrow=which(spec.dat$SpecID=="clav")
dat$AT_clav[which(dat$Te_clav >= spec.dat[wrow, "Tb20"] & dat$Te_clav <= spec.dat[wrow, "Tb80"])]=1
#pell
wrow=which(spec.dat$SpecID=="pell")
dat$AT_pell[which(dat$Te_pell >= spec.dat[wrow, "Tb20"] & dat$Te_pell <= spec.dat[wrow, "Tb80"])]=1
#dodg
wrow=which(spec.dat$SpecID=="dodg")
dat$AT_dodg[which(dat$Te_dodg >= spec.dat[wrow, "Tb20"] & dat$Te_dodg <= spec.dat[wrow, "Tb80"])]=1
#sang
wrow=which(spec.dat$SpecID=="sang")
dat$AT_sang[which(dat$Te_sang >= spec.dat[wrow, "Tb20"] & dat$Te_sang <= spec.dat[wrow, "Tb80"])]=1

#ALLOW TO USE SHADE
#clav
wrow=which(spec.dat$SpecID=="clav")
dat$AT_clav[which(dat$Te_clav >= spec.dat[wrow, "Tb20"] & dat$Te_clav <= spec.dat[wrow, "Tb80"])]=1
#pell
wrow=which(spec.dat$SpecID=="pell")
dat$AT_pell[which(dat$Te_pell >= spec.dat[wrow, "Tb20"] & dat$Te_pell <= spec.dat[wrow, "Tb80"])]=1
#dodg
wrow=which(spec.dat$SpecID=="dodg")
dat$AT_dodg[which(dat$Te_dodg >= spec.dat[wrow, "Tb20"] & dat$Te_dodg <= spec.dat[wrow, "Tb80"])]=1
#sang
wrow=which(spec.dat$SpecID=="sang")
dat$AT_sang[which(dat$Te_sang >= spec.dat[wrow, "Tb20"] & dat$Te_sang <= spec.dat[wrow, "Tb80"])]=1

#clav
wrow=which(spec.dat$SpecID=="clav")
dat$AT_clav[which(dat$Te_clav.shade >= spec.dat[wrow, "Tb20"] & dat$Te_clav <= spec.dat[wrow, "Tb80"])]=1
#pell
wrow=which(spec.dat$SpecID=="pell")
dat$AT_pell[which(dat$Te_pell.shade >= spec.dat[wrow, "Tb20"] & dat$Te_pell <= spec.dat[wrow, "Tb80"])]=1
#dodg
wrow=which(spec.dat$SpecID=="dodg")
dat$AT_dodg[which(dat$Te_dodg.shade >= spec.dat[wrow, "Tb20"] & dat$Te_dodg <= spec.dat[wrow, "Tb80"])]=1
#sang
wrow=which(spec.dat$SpecID=="sang")
dat$AT_sang[which(dat$Te_sang.shade >= spec.dat[wrow, "Tb20"] & dat$Te_sang <= spec.dat[wrow, "Tb80"])]=1

#clav
wrow=which(spec.dat$SpecID=="clav")
dat$AT_clav[which(dat$Te_clav.shade >= spec.dat[wrow, "Tb20"] & dat$Te_clav.shade <= spec.dat[wrow, "Tb80"])]=1
#pell
wrow=which(spec.dat$SpecID=="pell")
dat$AT_pell[which(dat$Te_pell.shade >= spec.dat[wrow, "Tb20"] & dat$Te_pell.shade <= spec.dat[wrow, "Tb80"])]=1
#dodg
wrow=which(spec.dat$SpecID=="dodg")
dat$AT_dodg[which(dat$Te_dodg.shade >= spec.dat[wrow, "Tb20"] & dat$Te_dodg.shade <= spec.dat[wrow, "Tb80"])]=1
#sang
wrow=which(spec.dat$SpecID=="sang")
dat$AT_sang[which(dat$Te_sang.shade >= spec.dat[wrow, "Tb20"] & dat$Te_sang.shade <= spec.dat[wrow, "Tb80"])]=1

#clav
wrow=which(spec.dat$SpecID=="clav")
dat$AT_clav[which(dat$Te_clav >= spec.dat[wrow, "Tb20"] & dat$Te_clav.shade <= spec.dat[wrow, "Tb80"])]=1
#pell
wrow=which(spec.dat$SpecID=="pell")
dat$AT_pell[which(dat$Te_pell >= spec.dat[wrow, "Tb20"] & dat$Te_pell.shade <= spec.dat[wrow, "Tb80"])]=1
#dodg
wrow=which(spec.dat$SpecID=="dodg")
dat$AT_dodg[which(dat$Te_dodg >= spec.dat[wrow, "Tb20"] & dat$Te_dodg.shade <= spec.dat[wrow, "Tb80"])]=1
#sang
wrow=which(spec.dat$SpecID=="sang")
dat$AT_sang[which(dat$Te_sang >= spec.dat[wrow, "Tb20"] & dat$Te_sang.shade <= spec.dat[wrow, "Tb80"])]=1

#----------------------------------------
#PLOTS
#FIG 5

#Compare Te, MR, AT
#Aggregate metrics
par(mfrow=c(3,1), cex=1.2, bty="l",lwd=2,tcl=0.5, mar=c(1.5, 3, 0.4, 0.2), mgp=c(1.5, 0.25, 0), oma=c(2,0,0,0))
#cols=c("orange", "green", "purple", "black")
cols=c("black","black","black", "grey")
ltys=c("solid", "dotted", "dashed", "solid")

#Subset to August
dat.dt<- strptime(dat$datetime, format="%m/%d/%Y %H:%M", tz="MST") #as.Date
dataug= subset(dat, dat.dt$mon==7)
 
#subset to day
dataug= subset(dataug, dataug$Rad>5)

dat3= aggregate(dataug, by=list(dataug$site, dataug$J), FUN=mean, na.rm=TRUE)
dat3.sd= aggregate(dataug, by=list(dataug$site, dataug$J), FUN=sd, na.rm=TRUE)
dat3.count= aggregate(dataug, by=list(dataug$site, dataug$J), FUN=count)
dat3.se=dat3.sd
dat3.se[,7:27] = dat3.sd[,7:27]/sqrt(dat3.count[,7:27])

dat2= aggregate(dataug, by=list(dataug$site), FUN=mean, na.rm=TRUE)
dat2.sd= aggregate(dataug, by=list(dataug$site), FUN=sd, na.rm=TRUE)
dat2.count= aggregate(dataug, by=list(dataug$site), FUN=count)
dat2.se=dat2.sd
dat2.se[,7:47] = dat2.sd[,7:47]/sqrt(dat2.count[,7:47])

dat2=dat2[order(dat2$Elev),]
dat2.se=dat2.se[order(dat2$Elev),]

#plot Te
plot(dat2$Elev, dat2$Te_clav, col=cols[1], ylim=range(33,43), xlim=range(1600,3100), type="b", xlab="", ylab="Te (?C)", lty=ltys[1])
points(dat2$Elev, dat2$Te_dodg, col=cols[2], type="b", lty=ltys[2])
points(dat2$Elev, dat2$Te_pell, col=cols[3], type="b", lty=ltys[3])
points(dat2$Elev, dat2$Te_sang, col=cols[4], type="b", lty=ltys[4])
arrows(dat2$Elev, dat2$Te_clav-dat2.se$Te_clav, dat2$Elev, dat2$Te_clav+dat2.se$Te_clav, code=3, angle=90,length=0.1, col=cols[1])
arrows(dat2$Elev, dat2$Te_dodg-dat2.se$Te_dodg, dat2$Elev, dat2$Te_dodg+dat2.se$Te_dodg, code=3, angle=90,length=0.1, col=cols[2])
arrows(dat2$Elev, dat2$Te_pell-dat2.se$Te_pell, dat2$Elev, dat2$Te_pell+dat2.se$Te_pell, code=3, angle=90,length=0.1, col=cols[3])
arrows(dat2$Elev, dat2$Te_sang-dat2.se$Te_sang, dat2$Elev, dat2$Te_sang+dat2.se$Te_sang, code=3, angle=90,length=0.1, col=cols[4])

spec.lab=c("A. clavatus", "M. dodgei", "C. pellucida", "M. sanguinipes")
legend("top", spec.lab, pch=1, col=cols,lty=ltys, bty="n", ncol=2, cex=0.8)

#PLOT Te SHADE
#plot(dat2$Elev, dat2$Te_clav.shade, col=cols[1], ylim=range(23,35), xlim=range(1600,3100), type="b", xlab="", ylab="Te (?C)")
#points(dat2$Elev, dat2$Te_dodg.shade, col=cols[2], type="b")
#points(dat2$Elev, dat2$Te_pell.shade, col=cols[3], type="b")
#points(dat2$Elev, dat2$Te_sang.shade, col=cols[4], type="b")

#-----------
##Plot Ta and radiation
#par(mfrow=c(2,1), cex=1.2, bty="l",lwd=2,tcl=0.5, mar=c(1.5, 3, 0.4, 0.2), mgp=c(1.5, 0.25, 0), oma=c(2,0,0,0))

##plot Ta
#plot(dat2$Elev, dat2$SoilTemp, col="black", ylim=range(25,35), xlim=range(1600,3100), type="b", xlab="Elevation (m)", ylab="Ta (?C)")
#arrows(dat2$Elev, dat2$SoilTemp-dat2.se$SoilTemp, dat2$Elev, dat2$SoilTemp+dat2.se$SoilTemp, code=3, angle=90,length=0.1, col=cols[4])

##Plot radiation
#plot(dat2$Elev, dat2$Rad, col="black", ylim=range(dat2$Rad), xlim=range(1600,3100), type="b", xlab="Elevation (m)", ylab="Radiation (W/m^2)")
##arrows(dat2$Elev, dat2$Rad-dat2.se$Rad, dat2$Elev, dat2$Rad+dat2.se$Rad, code=3, angle=90,length=0.1, col=cols[4])

#mtext("Elevation (m)", side = 1, line=2, cex=1.5)
#-----------

#plot MR
plot(dat2$Elev, dat2$MR_clav, col=cols[1], ylim=range(0,0.8), xlim=range(1600,3100), type="b", xlab="", ylab="MR (ml CO2/h)", lty=ltys[1])
points(dat2$Elev, dat2$MR_dodg, col=cols[2], type="b", lty=ltys[2])
points(dat2$Elev, dat2$MR_pell, col=cols[3], type="b", lty=ltys[3])
points(dat2$Elev, dat2$MR_sang, col=cols[4], type="b", lty=ltys[4])
arrows(dat2$Elev, dat2$MR_clav-dat2.se$MR_clav, dat2$Elev, dat2$MR_clav+dat2.se$MR_clav, code=3, angle=90,length=0.1, col=cols[1])
arrows(dat2$Elev, dat2$MR_dodg-dat2.se$MR_dodg, dat2$Elev, dat2$MR_dodg+dat2.se$MR_dodg, code=3, angle=90,length=0.1, col=cols[2])
arrows(dat2$Elev, dat2$MR_pell-dat2.se$MR_pell, dat2$Elev, dat2$MR_pell+dat2.se$MR_pell, code=3, angle=90,length=0.1, col=cols[3])
arrows(dat2$Elev, dat2$MR_sang-dat2.se$MR_sang, dat2$Elev, dat2$MR_sang+dat2.se$MR_sang, code=3, angle=90,length=0.1, col=cols[4])

#plot AT
plot(dat2$Elev, dat2$AT_clav, col=cols[1], ylim=range(0.25,0.75), xlim=range(1600,3100), type="b", xlab="", ylab="Proportion activity", lty=ltys[1])
points(dat2$Elev, dat2$AT_dodg, col=cols[2], type="b", lty=ltys[2])
points(dat2$Elev, dat2$AT_pell, col=cols[3], type="b", lty=ltys[3])
points(dat2$Elev, dat2$AT_sang, col=cols[4], type="b", lty=ltys[4])
arrows(dat2$Elev, dat2$AT_clav-dat2.se$AT_clav, dat2$Elev, dat2$AT_clav+dat2.se$AT_clav, code=3, angle=90,length=0.1, col=cols[1])
arrows(dat2$Elev, dat2$AT_dodg-dat2.se$AT_dodg, dat2$Elev, dat2$AT_dodg+dat2.se$AT_dodg, code=3, angle=90,length=0.1, col=cols[2])
arrows(dat2$Elev, dat2$AT_pell-dat2.se$AT_pell, dat2$Elev, dat2$AT_pell+dat2.se$AT_pell, code=3, angle=90,length=0.1, col=cols[3])
arrows(dat2$Elev, dat2$AT_sang-dat2.se$AT_sang, dat2$Elev, dat2$AT_sang+dat2.se$AT_sang, code=3, angle=90,length=0.1, col=cols[4])

mtext("Elevation (m)", side = 1, line=2, cex=1.5)

#----------------------------------------
#Sites in rows as a function of Date
par(mfrow=c(4,1), cex=1.2, lwd=2, mar=c(3, 3, 0.4, 0.2), mgp=c(1.5, 0.5, 0), oma=c(0,0,0,0))

#Subset to August
dat.dt<- strptime(dat$datetime, format="%m/%d/%Y %H:%M", tz="MST") #as.Date
dataug= subset(dat, dat.dt$mon==7)

#ELDORADO
dat1=subset(dataug, dataug$site=="Eldorado")
dat1$dt<- strptime(dat1$datetime, format="%m/%d/%Y %H:%M", tz="MST") #as.Date
plot(dat1$dt, dat1$Te_sang, col=cols[1], ylim=range(10,70), type="l", xlab="Date", ylab="Te (?C)")
points(dat1$dt, dat1$Te_dodg, col=cols[2], type="l")
points(dat1$dt, dat1$Te_pell, col=cols[3], type="l")
points(dat1$dt, dat1$Te_clav, col=cols[4], type="l")
  
dat1=subset(dataug, dataug$site=="A1")
dat1$dt<- strptime(dat1$datetime, format="%m/%d/%Y %H:%M", tz="MST") #as.Date
plot(dat1$dt, dat1$Te_sang, col=cols[1], ylim=range(10,70), type="l", xlab="Date", ylab="Te (?C)")
points(dat1$dt, dat1$Te_dodg, col=cols[2], type="l")
points(dat1$dt, dat1$Te_pell, col=cols[3], type="l")
points(dat1$dt, dat1$Te_clav, col=cols[4], type="l")

#B1
dat1=subset(dataug, dataug$site=="B1")
dat1$dt<- strptime(dat1$datetime, format="%m/%d/%Y %H:%M", tz="MST") #as.Date
plot(dat1$dt, dat1$Te_sang, col=cols[1], ylim=range(10,70), type="l", xlab="Date", ylab="Te (?C)")
points(dat1$dt, dat1$Te_dodg, col=cols[2], type="l")
points(dat1$dt, dat1$Te_pell, col=cols[3], type="l")
points(dat1$dt, dat1$Te_clav, col=cols[4], type="l")

#C1
dat1=subset(dataug, dataug$site=="C1")
dat1$dt<- strptime(dat1$datetime, format="%m/%d/%Y %H:%M", tz="MST") #as.Date
plot(dat1$dt, dat1$Te_sang, col=cols[1], ylim=range(10,70), type="l", xlab="Date", ylab="Te (?C)")
points(dat1$dt, dat1$Te_dodg, col=cols[2], type="l")
points(dat1$dt, dat1$Te_pell, col=cols[3], type="l")
points(dat1$dt, dat1$Te_clav, col=cols[4], type="l")

#----------------------------------------
par(mfcol=c(3,4), cex=1.2, lwd=2, mar=c(3, 3, 0.4, 0.2), mgp=c(1.5, 0.5, 0), oma=c(0,0,0,0))

for(i in 1:4){
if(i==1) dat1=subset(dat3, dat3$Group.1=="Eldorado")
if(i==2) dat1=subset(dat3, dat3$Group.1=="A1")
if(i==3) dat1=subset(dat3, dat3$Group.1=="B1")
if(i==4) dat1=subset(dat3, dat3$Group.1=="C1")

plot(dat1$J, dat1$Te_sang, col=cols[1], ylim=range(20,60), type="l", xlab="Date", ylab="Te (?C)")
points(dat1$J, dat1$Te_dodg, col=cols[2], type="l")
points(dat1$J, dat1$Te_pell, col=cols[3], type="l")
points(dat1$J, dat1$Te_clav, col=cols[4], type="l")

plot(dat1$J, dat1$MR_sang, col=cols[1], ylim=range(0,0.4), type="l", xlab="Date", ylab="MR")
points(dat1$J, dat1$MR_dodg, col=cols[2], type="l")
points(dat1$J, dat1$MR_pell, col=cols[3], type="l")
points(dat1$J, dat1$MR_clav, col=cols[4], type="l")

plot(dat1$J, dat1$AT_sang, col=cols[1], ylim=range(0,1), type="l", xlab="Date", ylab="AT")
points(dat1$J, dat1$AT_dodg, col=cols[2], type="l")
points(dat1$J, dat1$AT_pell, col=cols[3], type="l")
points(dat1$J, dat1$AT_clav, col=cols[4], type="l")
}
