################################################################################
##  codominance_analyses_figures.R: Merge datasets and perform analyses. Generate figures for publication.
##
##  Author: Kimberly Komatsu
##  Date created: June 9, 2020
################################################################################

library(tidyverse)
library(jmv)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\codominance\\data')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}


###read in data
GEx <- read.csv('GEx\\GEx_codominance_06112020.csv')%>%
  mutate(database=paste('GEx',trt))%>%
  mutate(method=ifelse(plot.size.shape %in% c('12.25','10 m long (50 points)','120 m long','1440 pin hits','15m long','20 m long','30 m line','60 m long','80 m long'), 'transect', 
                       ifelse(plot.size.shape %in% c('0.75m2','1 m2','1.25m2','1.4m2','1m2','2 m2','2.5m2','2-4m2','3.5 m2','4 m2','4m2','4m3','4m4','4m5','4m6','4m7','4m8','4m9','5 m2', '5m2'), 'small', 
                              ifelse(plot.size.shape %in% c('~23m2','10 m2','10m2','12m2','20m2','22.5m2','24m2','25m2','25m2?','8m2'), 'big', ifelse(is.na(plot.size.shape), 'NA', 'huge')))))

nutnet <- read.csv('nutnet\\NutNet_codominance_06092020.csv')%>%
  mutate(database='nutnet')

#nutnet evenness data
nutnetEven <- read.csv('nutnet\\nutnet_richEven_06122020.csv')%>%
  separate(exp_unit, into=c('site', 'block', 'plot', 'trt', 'year'), sep='::')%>%
  mutate(plot=as.integer(plot), year=as.integer(year), block=as.integer(block))%>%
  group_by(site, block, trt, year)%>%
  summarise(richness=mean(richness), Evar=mean(Evar))%>%
  ungroup()%>%
  group_by(site, trt, year)%>%
  summarise(richness=mean(richness), Evar=mean(Evar))%>%
  ungroup()%>%
  left_join(nutnet)

ggplot(data=subset(nutnetEven, year_trt==0), aes(x=Evar, y=num_codominants)) +
  geom_point() +
  xlab('Evar') + ylab('Number of Codominants')
  
ggplot(data=subset(nutnetEven, year_trt==0), aes(x=richness, y=num_codominants)) +
  geom_point() +
  xlab('Richness') + ylab('Number of Codominants')


#nutnet who are codom
nutnetFam <- read.csv('nutnet\\NutNet_codominants_list_06092020.csv')%>%
  filter(trt=='Control')%>%
  group_by(num_codominants, rank, Family)%>%
  summarise(count=length(Cmax))%>%
  ungroup()%>%
  mutate(family_group=ifelse(Family %in% c('Poaceae','Fabaceae','Compositae','Cyperaceae','Rubiaceae','Apiaceae','Geraniaceae'), as.character(Family), 'other'))

ggplot(nutnetFam, aes(x=rank, y=count, fill=family_group)) + 
  geom_bar(position="fill", stat="identity") +
  xlab('Rank') + ylab('Proportion') +
  facet_wrap(~num_codominants)

nutnetFunct <- read.csv('nutnet\\NutNet_codominants_list_06092020.csv')%>%
  filter(trt=='Control')%>%
  group_by(num_codominants, rank, functional_group)%>%
  summarise(count=length(Cmax))%>%
  ungroup()

ggplot(nutnetFunct, aes(x=rank, y=count, fill=functional_group)) + 
  geom_bar(position="fill", stat="identity") +
  xlab('Rank') + ylab('Proportion') +
  facet_wrap(~num_codominants)

# corre <- read.csv('CoRRE\\corre_codominance_06092020.csv')

###GEx methods
GEx$method <- factor(GEx$method, levels=c('transect','small','big','huge','NA'))
ggplot(data=barGraphStats(data=GEx, variable="num_codominants", byFactorNames=c("method",'database')), aes(x=method, y=mean)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  ylab('Number of Codominants') +
  facet_wrap(~database)
#export at 900x900


###NutNet map
library(wesanderson)
library("rnaturalearth")
library("rnaturalearthdata")
library("rgeos")
world <- ne_countries(scale = "medium", returnclass = "sf")

z.pal<-wes_palette("Zissou1", 100, type = "continuous")

ggplot(data=world)+
  theme(panel.background=element_rect(fill="aliceblue", color="aliceblue"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  geom_sf(color="black", fill="antiquewhite")+
  geom_point(data=nutnet, mapping=aes(x=longitude,y=latitude,fill=num_codominants),size=4.5,shape=21)+
  scale_fill_gradientn(colours = z.pal) +
  labs(fill="Number of Co-Dominants")+
  theme(legend.position = "top")+
  ylab(expression("Latitude "*degree*""))+
  xlab(expression("Longitude "*degree*""))
dev.off()


###NutNet families


###combine data
#only shared drivers
GExShared <- GEx%>%select(database, site, year, trt, Cmax, num_codominants, latitude, longitude, plant_gamma, MAT, temp_range, MAP, precip_cv, Ndep)
nutnetShared <- nutnet%>%select(database, site, year, trt, Cmax, num_codominants, latitude, longitude, plant_gamma, MAT, temp_range, MAP, precip_cv, Ndep)

codom <- rbind(GExShared,nutnetShared)%>%
  filter(!is.na(latitude), Ndep!='NULL')%>%
  mutate(Ndep=as.numeric(Ndep))

###how many sites have codominance?
ggplot(data=subset(codom, trt %in% c('G', 'U', 'Control')), aes(x=num_codominants)) +
  geom_histogram(color="black", fill="white", binwidth=1) +
  xlab('Number of Codominant Species') + ylab('Count') +
  facet_wrap(~database, scales='free_y')
#export at 2000x800



###MANCOVA
mancova(data=subset(codom, trt %in% c('G', 'U', 'Control') & num_codominants<15),
                                    deps=vars(Cmax, num_codominants),
                                    factors=c('database'),
                                    covs=c('plant_gamma', 'MAT', 'temp_range', 'MAP', 'precip_cv', 'Ndep'))

ggplot(data=subset(codom, trt %in% c('G', 'U', 'Control') & num_codominants<15), aes(x=plant_gamma, y=num_codominants, color=Cmax)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  xlab('Plant Gamma Diversity') + ylab('Number of Codominants') +
  facet_wrap(~database)
#export at 2000x800

ggplot(data=subset(codom, trt %in% c('G', 'U', 'Control') & num_codominants<15), aes(x=MAP, y=num_codominants, color=Cmax)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  xlab('Mean Annual Precipitation') + ylab('Number of Codominants') +
  facet_wrap(~database)
#export at 2000x800

#GEx only -- plot size or number?
ggplot(data=subset(GEx, trt %in% c('G', 'U') & num_codominants<15), aes(x=PlotSize, y=num_codominants, color=Cmax)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  xlab('Plot Size') + ylab('Number of Codominants') +
  facet_wrap(~database)

ggplot(data=subset(GEx, trt %in% c('G', 'U') & num_codominants<15), aes(x=Num_plots, y=num_codominants, color=Cmax)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  xlab('Plot Size') + ylab('Number of Codominants') +
  facet_wrap(~database)
#export at 2000x800









#code that needs modification:


###figures and models for number of codominants-------------
###multiple regression -- effects of grazing
#number of codominants - driven by species richness, biogeographic realm, marginally grazer presence
summary(lm(num_codominants ~ precip + sprich + biogeographic.realm + trt, data=CmaxDrivers))
#remove the outlier at 17 - driven by species richness, biogeographic realm
summary(lm(num_codominants ~ precip + sprich + biogeographic.realm + trt, data=subset(CmaxDrivers, num_codominants<17)))

#figures - number of codominants 
ggplot(data=CmaxDrivers, aes(x=sprich, y=num_codominants, color=Cmax)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  xlab('Gamma Diversity') + ylab('Number of Codominants')
#export at 1000x800

ggplot(data=barGraphStats(data=CmaxDrivers, variable="num_codominants", byFactorNames=c("biogeographic.realm")), aes(x=biogeographic.realm, y=mean)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Biogeographic Realm') + ylab('Number of Codominants') +
  annotate('text', x=1, y=2.7, label='a', size=12) +
  annotate('text', x=2, y=2.5, label='a', size=12) +
  annotate('text', x=3, y=2.2, label='ab', size=12) +
  annotate('text', x=4, y=2.2, label='ab', size=12) +
  annotate('text', x=5, y=2.0, label='ab', size=12) +
  annotate('text', x=6, y=1.7, label='b', size=12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#export at 1000x1000

ggplot(data=barGraphStats(data=CmaxDrivers, variable="num_codominants", byFactorNames=c("trt")), aes(x=trt, y=mean)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Grazer Presence') + ylab('Number of Codominants')
#export at 800x800


###multiple regression -- within grazed plots, effects of herbivores
#number of codominants - driven by species richness, biogeographic realm
summary(lm(num_codominants ~ precip + sprich + biogeographic.realm + grazing.pressure + herbivore.spp, data=subset(CmaxDrivers, trt=='G')))
#remove the outlier at 17 - driven by species richness, biogeographic realm, number of herbivore species
summary(lm(num_codominants ~ precip + sprich + biogeographic.realm + grazing.pressure +herbivore.spp, data=subset(CmaxDrivers, num_codominants<17 & trt=='G')))

#figures - number of codominants 
ggplot(data=subset(CmaxDrivers, trt=='G'), aes(x=sprich, y=num_codominants, color=Cmax)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  xlab('Gamma Diversity') + ylab('Number of Codominants')
#export at 1000x800

ggplot(data=barGraphStats(data=subset(CmaxDrivers, trt=='G'), variable="num_codominants", byFactorNames=c("biogeographic.realm")), aes(x=biogeographic.realm, y=mean)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Biogeographic Realm') + ylab('Number of Codominants') +
  annotate('text', x=1, y=3.5, label='a', size=12) +
  annotate('text', x=2, y=2.8, label='ab', size=12) +
  annotate('text', x=3, y=2.5, label='ab', size=12) +
  annotate('text', x=4, y=2.4, label='b', size=12) +
  annotate('text', x=5, y=2.2, label='b', size=12) +
  annotate('text', x=6, y=1.7, label='b', size=12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#export at 1000x1000

ggplot(data=subset(CmaxDrivers, trt=='G'), aes(x=herbivore.spp, y=num_codominants, color=Cmax)) +
  geom_point(size=5) +
  geom_smooth(data=subset(CmaxDrivers, trt=='G'&num_codominants<17), method='lm') +
  xlab('Number Herbivore Species') + ylab('Number of Codominants')
#export at 1000x800


###figures and models for Cmax-------------
###multiple regression -- effects of grazing
#number of codominants - driven by precip, species richness, marginally grazer presence
summary(lm(Cmax ~ precip + sprich + biogeographic.realm + trt, data=CmaxDrivers))

#figures - Cmax
ggplot(data=CmaxDrivers, aes(x=precip, y=Cmax, color=Cmax)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  xlab('Mean Annual Precipitation') + ylab('Cmax') +
  theme(legend.position='none')
#export at 1000x800

ggplot(data=CmaxDrivers, aes(x=sprich, y=Cmax, color=Cmax)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  xlab('Gamma Diversity') + ylab('Cmax') +
  theme(legend.position='none')
#export at 1000x800

ggplot(data=barGraphStats(data=CmaxDrivers, variable="Cmax", byFactorNames=c("trt")), aes(x=trt, y=mean)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Grazer Presence') + ylab('Cmax') +
  annotate('text', x=1, y=31, label='a', size=12) +
  annotate('text', x=2, y=34, label='b', size=12)
#export at 800x800


###multiple regression -- within grazed plots, effects of herbivores
#number of codominants - driven by precip, species richness, biogeographic realm
summary(lm(Cmax ~ precip + sprich + biogeographic.realm + grazing.pressure + herbivore.spp, data=subset(CmaxDrivers, trt=='G')))


#figures - Cmax
ggplot(data=subset(CmaxDrivers, trt=='G'), aes(x=sprich, y=Cmax, color=Cmax)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  xlab('Gamma Diversity') + ylab('Cmax') +
  theme(legend.position='none')
#export at 1000x800

ggplot(data=subset(CmaxDrivers, trt=='G'), aes(x=sprich, y=Cmax, color=Cmax)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  xlab('Mean Annual Precipitation') + ylab('Cmax') +
  theme(legend.position='none')
#export at 1000x800

ggplot(data=barGraphStats(data=subset(CmaxDrivers, trt=='G'), variable="Cmax", byFactorNames=c("biogeographic.realm")), aes(x=biogeographic.realm, y=mean)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Biogeographic Realm') + ylab('Cmax') +
  annotate('text', x=1, y=31, label='a', size=12) +
  annotate('text', x=2, y=29, label='a', size=12) +
  annotate('text', x=3, y=41, label='ab', size=12) +
  annotate('text', x=4, y=30, label='ab', size=12) +
  annotate('text', x=5, y=37, label='ab', size=12) +
  annotate('text', x=6, y=39, label='b', size=12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#export at 1000x1000
