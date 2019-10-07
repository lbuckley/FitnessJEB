#Explore survival and fecundity data in meta-analyses
library(ggplot2)
library(reshape2)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/meta_data")

#===============================================
#Anderson
#https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.13693
#binary fitness components (germination, seedling establishment, juvenile survival and adult survival
#Fecundity data were available at the level of the individual plant (e.g. number of flowers, fruits or seeds per individual) for 55 records (species by study combinations). For the remaining 40 records, fecundity data were presented on a population level, either per unit area (e.g. number of seeds or seed biomass per m2) or per plot (e.g. fecundity per experimental block). 

#some matches between s1 (survival) and s4 (fecundity)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/meta_data/Anderson_PlantMeta/")
s1= read.csv("S1.csv")
s4= read.csv("S4.csv")

#===========================================
#LEE-YAW
# no values

#===============================================
#Hargraves

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/meta_data/Hargreaves/")
fit= read.csv("Hargreaves_metaperfdata.csv")

#parameter
#Em=emergence (germination to young life stage for plants, hatching for insects), G=growth, S=survival, R=reproduction, 
#Ov=overall fitness=combination of >1 fitness component (Em, S, R) or lifetime fitness as presented in article (see Hargreaves et al 2013 for full explanation)

#RESTRICT to S, R, Ov
fit= fit[fit$parameter %in% c("S","R"),]
fit= fit[fit$limit %in% c("high","low"),] #"equitoral","polar",

#limit type:	divides range limits into G=geographic, E=elevational
#limit:	location of RL relative to the species entire range. high=high elev RL, low=low elev RL, polar=high latitude RL (northern RL in northern hemisphere, southern RL in southern hemisphere), equatorial=low latitude RL, east=eastern RL (longitudinal), west=western RL (longitudinal).

#melt
f.long= melt(data = fit, id.vars = c("reference","species","organism","limit.type","limit","parameter"), measure.vars = c("i.I", "i.E", "i.B", "e.I", "e.E", "e.B"))
f.long$home<-"home"
f.long$home[f.long$variable %in% c("i.E","e.I" )] <-"away"
f.long$home[f.long$variable %in% c("i.B","e.B" )] <-"beyond"
f.long$source<-"int"
f.long$source[f.long$variable %in% c("e.I","e.E","e.B")] <-"edge"
f.long$site<-"int"
f.long$site[f.long$variable %in% c("i.E","e.E")] <-"edge"
f.long$site[f.long$variable %in% c("i.B","e.B")] <-"beyond"
f.long$site= factor(f.long$site, levels=c("int","edge","beyond"))
f.long$group= paste(f.long$reference, f.long$species, f.long$limit, sep="_")

#----
#PLOT
ggplot(data=f.long[f.long$parameter %in% c("R","G"),], aes(x=site, y = value, color=source, group=group))+geom_point()+facet_grid(parameter~limit, scales="free_y") +
  stat_summary(fun.y=mean, geom="line")
ggplot(data=f.long[f.long$parameter=="S",], aes(x=site, y = value, color=source, group=group))+geom_point()+facet_grid(parameter~limit, scales="free_y") +
  stat_summary(fun.y=mean, geom="line")

ggplot(data=f.long[f.long$parameter %in% c("R","G"),], aes(x=site, y = log(value), color=source))+geom_boxplot()+facet_grid(~limit, scales="free_y")
ggplot(data=f.long[f.long$parameter %in% c("S"),], aes(x=site, y = value, color=source))+geom_boxplot()+facet_grid(~limit, scales="free_y")

#-----------------
#performance difference away to home
ggplot(data=fit, aes(x=parameter, y = i.EI, color=limit))+geom_boxplot()+facet_wrap(~reference, ncol=4) #, scales="free_y"
ggplot(data=fit, aes(x=parameter, y = e.IE, color=limit))+geom_boxplot()+facet_wrap(~reference, ncol=4) #, scales="free_y"


ggplot(data=fit, aes(x=parameter, y = i.BI, color=organism, group=parameter))+geom_boxplot()+facet_wrap(~limit, ncol=1) #, scales="free_y"
ggplot(data=fit, aes(x=parameter, y = i.EI, color=organism, group=parameter))+geom_boxplot()+facet_wrap(~limit, ncol=1) #, scales="free_y"





