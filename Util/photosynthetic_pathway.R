################################################################################
##  merge_community_diff_metrics.R: Mergine files of community difference metrics.
##
##  Author: Kimberly Komatsu and Sally Koerner
##  Date created: April 25, 2019
################################################################################
##Update by S Koerner on June 12, 2020
################################################################################

library(tidyverse)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\data')
#Sally's laptop
setwd("~/Dropbox/GEx_VirtualWorkshop_June2020")

###read in all data
familyList <- read.csv('GEx_sppfam_FINAL_10June2020.csv')%>% #list of all species and their families/photosynthetic pathways
  select(-genus_species, -tnrs_accepted_name)%>%
  filter(!is.na(clean_ejf)) %>% 
  rename(genus_species_clean=clean_ejf)

familyList2<-familyList %>% 
  group_by(genus_species_clean, family_ejf, pathway) %>% 
  summarise(num=n()) %>% 
  select(-num) %>% 
  ungroup() 

coverData <- read.csv('GEx_cleaned_11June2020.csv') %>% 
  mutate(drop=ifelse(genus_species_use=="#N/A"|genus_species_use=="Dead unidentified"|genus_species_use=="Leaf.Litter"|genus_species_use=="cactus__dead_", 1, 0))%>%
  filter(drop!=1) %>% 
  select(-drop)%>% 
  left_join(familyList2)%>%
  unique()

photopath <- coverData%>%
  mutate(pathway=ifelse(is.na(pathway), "unknown", as.character(pathway)))%>%
  group_by(site, block, plot, trt, pathway)%>%
  summarise(photo_path=sum(relcov))%>%
  ungroup()%>%
  spread(key=pathway, value=photo_path, fill=0)%>%
  gather(key=pathway, value=relcov, C3:unknown)

#calculate at the site level
photopathSite <- photopath%>%
  group_by(site, block, pathway)%>%
  summarise(relcov=mean(relcov))%>%
  ungroup()%>%
  group_by(site, pathway)%>%
  summarise(relcov=mean(relcov))%>%
  ungroup()%>%
  spread(key=pathway, value=relcov, fill=0)


write.csv(photopathSite, 'percent_photosynthetic_pathway.csv', row.names=F)
