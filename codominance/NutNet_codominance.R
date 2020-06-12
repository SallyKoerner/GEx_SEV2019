################################################################################
##  NutNet_codominance.R: Calculate codominance of a community for NutNet database.
##
##  Author: Kimberly Komatsu
##  Date created: June 9, 2020
################################################################################

library(psych)
library(codyn)
library(tidyverse)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\codominance\\data\\nutnet')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###read in data
nutnet <- read.csv('full-cover-03-June-2020.csv')%>%
  rename(cover=max_cover, genus_species=Taxon)%>%
  mutate(exp_unit=paste(site_code, block, plot, trt, year, sep='::'))%>%
  filter(!(genus_species %in% c('GROUND', 'OTHER LITTER', 'OTHER ARISTIDA CONTORTA (DEAD)', 'OTHER SALSOLA KALI (DEAD)', 'OTHER TRIODIA BASEDOWII (DEAD)')))

#############################################
#####calculate Cmax (codominance metric)#####

#calculate relative abundance
relCover <- nutnet%>%
  group_by(exp_unit)%>%
  summarise(totcov=sum(cover))%>%
  ungroup()%>%
  right_join(nutnet)%>%
  mutate(relcov=(cover/totcov)*100)%>%
  select(-cover, -totcov)

evenness <- relCover%>%
  community_structure(time.var = 'year', abundance.var = 'relcov',
                      replicate.var = 'exp_unit', metric = c("Evar", "SimpsonEvenness", "EQ"))

# write.csv(evenness, 'nutnet_richEven_06122020.csv')

#generate rank of each species in each plot by relative cover, with rank 1 being most abundant
rankOrder <- relCover%>%
  group_by(exp_unit)%>%
  mutate(rank = as.numeric(order(order(relcov, decreasing=TRUE))))%>%
  ungroup()

###calculating harmonic means for all subsets of rank orders
#make a new dataframe with just the label
expUnit=nutnet%>%
  select(exp_unit)%>%
  unique()

#makes an empty dataframe
harmonicMean=data.frame(row.names=1) 

### NOTE: this code takes about 30 mins to run, so use the output in the dropbox unless there is a reason to re-run it
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
  mutate(plot=as.integer(plot), year=as.integer(year), block=as.integer(block))

codomSppList <- Cmax%>%
  left_join(rankOrder)%>%
  group_by(exp_unit)%>%
  filter(rank<=num_codominants)%>%
  ungroup()

#write.csv(codomSppList, 'NutNet_codominants_list_06092020.csv', row.names=F)

#histogram of codom
ggplot(data=subset(codomSppList, year_trt==0), aes(x=num_codominants)) +
  geom_histogram(color="black", fill="white", binwidth=1) +
  xlab('Number of Codominant Species') + ylab('Count')
#export at 1000x800

ggplot(data=subset(codomSppList, year_trt==0), aes(x=Cmax, y=num_codominants)) +
  geom_point() +
  xlab('Cmax') + ylab('Number of Codominants')
#export at 800x800


###what drives codominance?

#read in site-level data
siteData <- read.csv('comb-by-plot-clim-soil-diversity-03-Jun-2020.csv')%>%
  select(site_code, continent, country, region, managed, burned, grazed, anthropogenic, habitat, elevation, latitude, longitude, site_richness, MAT_v2, ANN_TEMP_RANGE_v2, MAP_v2, MAP_VAR_v2, N_Dep, experiment_type, year, year_trt, first_nutrient_year, first_fenced_year, site_native_richness, site_introduced_richness)%>%
  unique()%>%
  rename(site=site_code, plant_gamma=site_richness, MAT=MAT_v2, temp_range=ANN_TEMP_RANGE_v2, MAP=MAP_v2, precip_cv=MAP_VAR_v2, Ndep=N_Dep)

#get site-level average cmax and number of codominants
CmaxDrivers <- Cmax%>%
  group_by(site, year, trt, block)%>%
  summarise(num_codominants=mean(num_codominants), Cmax=mean(Cmax))%>%
  ungroup()%>%
  group_by(site, year, trt)%>%
  summarise(num_codominants=mean(num_codominants), Cmax=mean(Cmax))%>%
  ungroup()%>%
  left_join(siteData)

# write.csv(CmaxDrivers, 'NutNet_codominance_06092020.csv', row.names=F)

ggplot(data=CmaxDrivers, aes(x=Cmax, y=num_codominants)) +
  geom_point() +
  xlab('Cmax') + ylab('Number of Codominants')
#export at 800x800
