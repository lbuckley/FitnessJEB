#Explore survival and fecundity data in meta-analyses
library(ggplot2)
library(reshape2)
library(Hmisc)
library(cowplot)

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
#restrict to warming treatments
s1= s1[s1$treatment %in% c("infrared_heater","warming"),] #"snow_removal",
s4= s4[s4$treatment %in% c("active_warming","passive_warming"),]

s1$surv_e= s1$ai/s1$n1i
s1$surv_c= s1$ci/s1$n2i
s1$surv_d= (s1$surv_e-s1$surv_c)/s1$surv_c

s4$fit_d= (s4$m1i-s4$m2i)/s4$m2i

#survival
s.plot= ggplot(data=s1, aes(x=log(elevation.m.), y = surv_d, group=1, color=treatment))+geom_point()+ylim(-1,2)+
  geom_smooth(method="lm",color="black",se=FALSE)+theme(legend.position = "bottom")+geom_hline(yintercept=0)+
  ylab("proportional survival change") +xlab("log(elevation) (m)")
#omits some data

#fedundity
f.plot=ggplot(data=s4, aes(x=log(elevation.m.), y = fit_d, group=1, color=treatment))+geom_point()+ylim(-1,2)+
  geom_smooth(method="lm",color="black",se=FALSE)+theme(legend.position = "bottom")+geom_hline(yintercept=0)+
  ylab("proportional fecundity change") +xlab("log(elevation) (m)")  
#omits some data

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/figures/meta_data/")
pdf("AndersonFig.pdf", height = 8, width = 10)
plot_grid(s.plot,f.plot)
dev.off()

#===============================================
#Hargraves

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/meta_data/Hargreaves/")
fit= read.csv("Hargreaves_metaperfdata.csv")

#parameter
#Em=emergence (germination to young life stage for plants, hatching for insects), G=growth, S=survival, R=reproduction, 
#Ov=overall fitness=combination of >1 fitness component (Em, S, R) or lifetime fitness as presented in article (see Hargreaves et al 2013 for full explanation)

#RESTRICT to S, R, Ov
fit= fit[fit$parameter %in% c("S","R","G"),]
fit= fit[fit$limit %in% c("high","low","equitoral","polar"),] 

#limit type:	divides range limits into G=geographic, E=elevational
#limit:	location of RL relative to the species entire range. high=high elev RL, low=low elev RL, polar=high latitude RL (northern RL in northern hemisphere, southern RL in southern hemisphere), equatorial=low latitude RL, east=eastern RL (longitudinal), west=western RL (longitudinal).

# #UNNORMALIZED
# #melt
# f.long= melt(data = fit, id.vars = c("reference","species","organism","limit.type","limit","parameter"), measure.vars = c("i.I", "i.E", "i.B", "e.I", "e.E", "e.B"))
# f.long$home<-"home"
# f.long$home[f.long$variable %in% c("i.E","e.I" )] <-"away"
# f.long$home[f.long$variable %in% c("i.B","e.B" )] <-"beyond"
# f.long$source<-"int"
# f.long$source[f.long$variable %in% c("e.I","e.E","e.B")] <-"edge"
# f.long$site<-"int"
# f.long$site[f.long$variable %in% c("i.E","e.E")] <-"edge"
# f.long$site[f.long$variable %in% c("i.B","e.B")] <-"beyond"

#-------------
#NORMALIZED
#melt
f.long= melt(data = fit, id.vars = c("reference","species","organism","limit.type","limit","parameter"), measure.vars = c("i.BI", "i.EI", "i.BE", "e.BI", "e.IE", "e.BE"))
#drop extra comparisons
f.long= f.long[f.long$variable %in% c("i.EI","i.BI","e.IE"), ] #drop edge to beyond: "e.BE"

f.long$source<-"interior"
f.long$source[f.long$variable %in% c("e.IE","e.BE")] <-"edge"
f.long$site<-NA
f.long$site[f.long$variable %in% c("e.IE")] <-"interior"
f.long$site[f.long$variable %in% c("i.EI")] <-"edge"
f.long$site[f.long$variable %in% c("i.BI", "e.BE")] <-"beyond"

#other variables
f.long$site= factor(f.long$site, levels=c("int","edge","beyond"))
f.long$group= paste(f.long$reference, f.long$species, f.long$limit, sep="_")
f.long$edge<-"cold limit"
f.long$edge[f.long$limit %in% c("low","equitoral")] <-"warm limit"
#change parameter name
f.long$parameter.lab<-NA
f.long$parameter.lab[f.long$parameter=="G"]<-"growth"
f.long$parameter.lab[f.long$parameter=="R"]<-"reproduction"
f.long$parameter.lab[f.long$parameter=="S"]<-"survival"

#---------------------------
#FIGURE

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/figures/meta_data/")
pdf("HargreavesFig.pdf", height = 8, width = 10)
ggplot(data=f.long, aes(x=site, y = value, color=parameter.lab, shape=source))+geom_jitter(size=2, width=0.3, height=0)+facet_grid(edge~parameter.lab) +
  geom_hline(yintercept=0)+stat_summary(fun.data=mean_sdl,fun.args = list(mult=1), geom="errorbar",width=0.2, color="black")+
  stat_summary(fun.y=mean, geom="point", color="black", size=4, shape=0) +theme(legend.position = "bottom")+guides(colour=FALSE)+ #+theme(strip.text.y = element_text(angle = 360))+stat_summary(fun.y=mean, geom="point")
  ylim(-2.5,2.5)
  dev.off()

#==============================
#ANIMALS
  
#from Hargreaves
#check animals
  ggplot(data=f.long[f.long$organism=="Animal",], aes(x=site, y = value, color=parameter.lab, shape=source))+geom_jitter(size=2, width=0.3, height=0)+facet_grid(edge+reference~parameter.lab) +
    geom_hline(yintercept=0)+stat_summary(fun.data=mean_sdl,fun.args = list(mult=1), geom="errorbar",width=0.2, color="black")+
    stat_summary(fun.y=mean, geom="point", color="black", size=4, shape=0) +theme(legend.position = "bottom")+guides(colour=FALSE)+ #+theme(strip.text.y = element_text(angle = 360))+stat_summary(fun.y=mean, geom="point")
    ylim(-2.5,2.5)
  
  
  


