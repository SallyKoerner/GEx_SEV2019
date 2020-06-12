################################################################################
##  GEx_codominance.R: Calculate codominance of a community for GEx database.
##
##  Author: Kimberly Komatsu
##  Date created: April 24, 2019
################################################################################

library(psych)
library(tidyverse)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\codominance\\data\\GEx')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###read in data
GEx <- read.csv('GEx_cleaned_11June2020.csv')%>%
  mutate(drop=ifelse(genus_species_use=="#N/A"|genus_species_use=="Dead unidentified"|genus_species_use=="Leaf.Litter"|genus_species_use=="cactus__dead_", 1, 0))%>%
  filter(drop!=1)%>%
  select(-genus_species, -genus_species_clean, -drop)%>%
  rename(genus_species=genus_species_use, cover=relcov)%>%
  group_by(site, year, exage, block, trt, genus_species)%>%
  summarise(cover=mean(cover))%>%
  ungroup()%>%
  mutate(exp_unit=paste(site, block, trt, year, sep='::'))

  
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
  separate(exp_unit2, into=c('site', 'block', 'trt', 'year'), sep='::')%>%
  mutate(year=as.integer(year))

codomSppList <- Cmax%>%
  left_join(rankOrder)%>%
  group_by(exp_unit)%>%
  filter(rank<=num_codominants)%>%
  ungroup()

# write.csv(codomSppList, 'GEx_codominants_list_06112020.csv', row.names=F)

#histogram of codom
ggplot(data=codomSppList, aes(x=num_codominants)) +
  geom_histogram(color="black", fill="white", binwidth=1) +
  xlab('Number of Codominant Species') + ylab('Count')
#export at 1000x800

ggplot(data=codomSppList, aes(x=Cmax, y=num_codominants)) +
  geom_point() +
  xlab('Cmax') + ylab('Number of Codominants')
#export at 800x800


###what drives codominance?

#read in site-level data
siteData <- read.csv('GEx-metadata-with-other-env-layers-v2.csv')%>%
  select(-X)%>%
  rename(plant_gamma=sprich, Ndep=N.deposition1993, latitude=Final.Lat, longitude=Final.Long, MAT=bio1, temp_range=bio7, MAP=bio12, precip_cv=bio15)

#get site-level average cmax and number of codominants
CmaxDrivers <- Cmax%>%
  group_by(site, year, trt)%>%
  summarise(num_codominants=mean(num_codominants), Cmax=mean(Cmax))%>%
  ungroup()%>%
  left_join(siteData)

# write.csv(CmaxDrivers, 'GEx_codominance_06112020.csv', row.names=F)

ggplot(data=CmaxDrivers, aes(x=Cmax, y=num_codominants)) +
  geom_point() +
  xlab('Cmax') + ylab('Number of Codominants')
#export at 800x800