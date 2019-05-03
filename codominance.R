################################################################################
##  codominance.R: Function to calculate codominance of a community.
##
##  Author: Kimberly Komatsu
##  Date created: April 24, 2019
################################################################################

library(psych)
library(tidyverse)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\data')


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
GEx <- read.csv('All_Cleaned_April2019_V2.csv')%>%
  rename(cover=relcov)%>%
  mutate(exp_unit=paste(site, block, plot, trt, year, sep='::'))

# #subset a test dataset of one plot to create a function
# test <- GEx%>%
#   filter(site=='Konza'&block=='A'&plot==1)


#############################################
#####calculate Cmax (codominance metric)#####

#calculate relative abundance
relCover <- GEx%>%
  group_by(exp_unit)%>%
  summarise(totcov=sum(cover))%>%
  ungroup()%>%
  right_join(GEx)%>%
  mutate(relcov=(cover/totcov)*100)%>%
  select(-cover, -totcov)

#generate rank of each species in each plot by relative cover, with rank 1 being most abundant
rankOrder <- relCover%>%
  group_by(exp_unit)%>%
  mutate(rank = as.numeric(order(order(relcov, decreasing=TRUE))))%>%
  ungroup()

###calculating harmonic means for all subsets of rank orders
#make a new dataframe with just the label
expUnit=GEx%>%
  select(exp_unit)%>%
  unique()

#makes an empty dataframe
harmonicMean=data.frame(row.names=1) 

#calculate harmonic means
for(i in 1:length(expUnit$exp_unit)) {
  
  #creates a dataset for each unique experimental unit
  subset <- rankOrder[rankOrder$exp_unit==as.character(expUnit$exp_unit[i]),]%>%
    select(exp_unit, genus_species, relcov, rank)
  
  for(j in 1:length(subset$rank)) {
    
    #creates a dataset for each series of ranks from 1 through end of the number of ranks
    subset2 <- subset[subset$rank<=j,]
    
    #calculate harmonic mean of values
    mean <- harmonic.mean(subset2$relcov)
    meanData <- data.frame(exp_unit=unique(subset2$exp_unit),
                           num_ranks=j, 
                           harmonic_mean=mean)
    
    harmonicMean=rbind(meanData, harmonicMean)
    
  }
  
}

differenceData <- harmonicMean%>%
  left_join(rankOrder)%>%
  filter(rank==num_ranks+1)%>% #only keep the next most abundant species after the number that went into the calculation of the harmonic mean
  mutate(difference=harmonic_mean-relcov) #calculates difference between harmonic mean and the relative cover of the next most abundant species
  
Cmax <- differenceData%>%
  group_by(exp_unit)%>%
  summarise(Cmax=max(difference))%>%
  ungroup()%>%
  left_join(differenceData)%>%
  filter(Cmax==difference)%>%
  rename(num_codominants=num_ranks)%>%
  select(exp_unit, Cmax, num_codominants)%>%
  mutate(exp_unit2=exp_unit)%>%
  separate(exp_unit2, into=c('site', 'block', 'plot', 'trt', 'year'), sep='::')%>%
  mutate(plot=as.integer(plot), year=as.integer(year))

codomSppList <- Cmax%>%
  left_join(rankOrder)%>%
  group_by(exp_unit)%>%
  filter(rank<=num_codominants)%>%
  ungroup()

#write.csv(codomSppList, 'GEx_codominane_04242019.csv', row.names=F)

#histogram of codom
ggplot(data=codomSppList, aes(x=num_codominants)) +
  geom_histogram(color="black", fill="white", binwidth=1) +
  xlab('Number of Codominant Species') + ylab('Count')


###what drives codominance?

#read in site-level data
siteData <- read.csv('Meta_SEV2019_v2.csv')

#get site-level average cmax and number of codominants
CmaxDrivers <- Cmax%>%
  group_by(site, block, trt)%>%
  summarise(num_codominants=mean(num_codominants), Cmax=mean(Cmax))%>%
  ungroup()%>%
  group_by(site, trt)%>%
  summarise(num_codominants=mean(num_codominants), Cmax=mean(Cmax))%>%
  ungroup()%>%
  left_join(siteData)


ggplot(data=CmaxDrivers, aes(x=Cmax, y=num_codominants)) +
  geom_point() +
  xlab('Cmax') + ylab('Number of Codominants')


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
#remove the outlier at 17 - driven by species richness, biogeographic realm, number of herbivore species
summary(lm(num_codominants ~ precip + sprich + biogeographic.realm + grazing.pressure +herbivore.spp, data=subset(CmaxDrivers, num_codominants<17 & trt=='G')))

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
