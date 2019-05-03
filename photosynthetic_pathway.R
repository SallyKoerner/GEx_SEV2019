################################################################################
##  merge_community_diff_metrics.R: Mergine files of community difference metrics.
##
##  Author: Kimberly Komatsu and Sally Koerner
##  Date created: April 25, 2019
################################################################################

library(tidyverse)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\data')

###read in all data
familyList <- read.csv('GEx_species_family_final_jn.csv')%>% #list of all species and their families/photosynthetic pathways
  select(-X, -genus_species, -tnrs_accepted_name)%>%
  filter(!is.na(clean_ejf))
coverData <- read.csv('GEx_cleaned_v3.csv')%>% #all cover data
  left_join(familyList)%>%
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


# write.csv(photopathSite, 'percent_photosynthetic_pathway.csv', row.names=F)
