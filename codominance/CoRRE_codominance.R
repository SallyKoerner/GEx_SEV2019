################################################################################
##  CoRRE_codominance.R: Calculate codominance of a community for CoRRE database.
##
##  Author: Kimberly Komatsu
##  Date created: June 12, 2020
################################################################################

library(psych)
library(tidyverse)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\codominance\\data\\CoRRE')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###read in data
corre <- read.csv('SpeciesRawAbundance_Nov2019.csv')%>%
  select(-X)%>%
  left_join(read.csv('corre2trykey.csv'))%>%
  select(-genus_species)%>%
  rename(genus_species=species_matched, cover=abundance, plot=plot_id, trt=treatment, year=calendar_year)%>%
  mutate(exp_unit=paste(site_code, project_name, plot, trt, year, sep='::'))

#############################################
#####calculate Cmax (codominance metric)#####

#calculate relative abundance
relCover <- corre%>%
  group_by(exp_unit)%>%
  summarise(totcov=sum(cover))%>%
  ungroup()%>%
  right_join(corre)%>%
  mutate(relcov=(cover/totcov)*100)%>%
  select(-cover, -totcov)

#generate rank of each species in each plot by relative cover, with rank 1 being most abundant
rankOrder <- relCover%>%
  group_by(exp_unit)%>%
  mutate(rank = as.numeric(order(order(relcov, decreasing=TRUE))))%>%
  ungroup()

###calculating harmonic means for all subsets of rank orders
#make a new dataframe with just the label
expUnit=corre%>%
  select(exp_unit)%>%
  unique()

#makes an empty dataframe
harmonicMean=data.frame(row.names=1) 

### NOTE: this code takes about 5 hours to run, so use the output in the dropbox unless there is a reason to re-run it
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
  separate(exp_unit2, into=c('site', 'project_name', 'plot', 'trt', 'year'), sep='::')%>%
  mutate(year=as.integer(year))

codomSppList <- Cmax%>%
  left_join(rankOrder)%>%
  group_by(exp_unit)%>%
  filter(rank<=num_codominants)%>%
  ungroup()

#write.csv(codomSppList, 'corre_codominants_list_08242020.csv', row.names=F)

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

#fix starting here...

#read in site-level data
siteData <- read.csv('SiteExperimentDetails_March2019.csv')%>%
  left_join(read.csv('ExperimentInformation_March2019.csv'))%>%
  select(-X)%>%
  # select(site_code, continent, country, region, managed, burned, grazed, anthropogenic, habitat, elevation, latitude, longitude, site_richness, MAT_v2, ANN_TEMP_RANGE_v2, MAP_v2, MAP_VAR_v2, N_Dep, experiment_type, year, year_trt, first_nutrient_year, first_fenced_year, site_native_richness, site_introduced_richness)%>%
  unique()%>%
  rename(site=site_code, plant_gamma=rrich)

#get site-level average cmax and number of codominants
CmaxDrivers <- Cmax%>%
  group_by(site, year, trt)%>%
  summarise(num_codominants=mean(num_codominants), Cmax=mean(Cmax))%>%
  ungroup()%>%
  left_join(siteData)

# write.csv(CmaxDrivers, 'corre_codominance_08242020.csv', row.names=F)

ggplot(data=CmaxDrivers, aes(x=Cmax, y=num_codominants)) +
  geom_point() +
  xlab('Cmax') + ylab('Number of Codominants')
#export at 800x800
