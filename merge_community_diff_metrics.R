################################################################################
##  merge_community_diff_metrics.R: Mergine files of community difference metrics.
##
##  Author: Kimberly Komatsu and Sally Koerner
##  Date created: April 25, 2019
################################################################################

library(tidyverse)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\data')

###import data
domDiff <- read.csv('OnlineVersion_withSimpD_Apr2018.csv')%>%select(-X) #dominance difference, with a bunch of site level covariates
compDiff <- read.csv('gex_multdiff_ave.csv') #composition difference
RACdiff <- read.csv('gex_RACdiff_ave_lnRR.csv') #RAC differences
domIDdiff <- read.csv('Dom_SameDiff_bySite_SEV_April2019.csv')%>%select(-X)%>%rename(IDdiff=num) #difference in identity of dominant spp
codomDiff <- read.csv('GEx_codominane_04242019.csv')%>% #codominance imported, need to calculate diff
  select(site, block, plot, trt, num_codominants)%>%
  unique()%>%
  spread(key=trt, value=num_codominants)%>%
  mutate(num_codom_diff=log(G/U))%>%
  filter(!is.na(num_codom_diff))%>%
  group_by(site, block)%>%
  summarise(num_codom_diff=mean(num_codom_diff))%>%
  ungroup()%>%
  group_by(site)%>%
  summarise(num_codom_diff=mean(num_codom_diff))%>%
  ungroup()
CmaxDiff <- read.csv('GEx_codominane_04242019.csv')%>% #Cmax imported, need to calculate diff
  select(site, block, plot, trt, Cmax)%>%
  unique()%>%
  spread(key=trt, value=Cmax)%>%
  mutate(Cmax_diff=log(G/U))%>%
  filter(!is.na(Cmax_diff))%>%
  group_by(site, block)%>%
  summarise(Cmax_diff=mean(Cmax_diff))%>%
  ungroup()%>%
  group_by(site)%>%
  summarise(Cmax_diff=mean(Cmax_diff))%>%
  ungroup()
photopath <- read.csv('percent_photosynthetic_pathway.csv') #photosynthetic pathways
phylorealms <- read.csv('PhyloRealms.csv')%>%select(-X) #phylogenetic realms
domRankDiff <- read.csv('sitelevel_domspecies_rankchange.csv')%>%select(-X) #rank change of dominant species
domPercentDiff <- read.csv('DomIdentityNumChangeLRR_bySite_SEV_April2019.csv') #percent abundance difference of dominant species
blockNum <- read.csv('GEx_cleaned_v3.csv')%>%
  select(site, block)%>%
  unique()%>%
  group_by(site)%>%
  summarise(num_blocks=length(block))%>%
  ungroup()
  

###merge together
compDiffSite <- domDiff%>%
  left_join(compDiff)%>%
  left_join(RACdiff)%>%
  left_join(domIDdiff)%>%
  left_join(CmaxDiff)%>%
  left_join(codomDiff)%>%
  left_join(photopath)%>%
  left_join(phylorealms)%>%
  left_join(domRankDiff)%>%
  left_join(domPercentDiff)%>%
  left_join(blockNum)%>%
  mutate(site_dom=(GDom+UDom)/2) #calculate average dominance between grazed and ungrazed areas

# write.csv(compDiffSite, 'community_difference_allmetrics_siteavg_04262019b.csv', row.names=F)


####Sallys Add in strats here --- reimport above write

###import data
All <- read.csv('community_difference_allmetrics_siteavg_04262019b.csv')
Climate <- read.csv('climatePCAs.csv')
Meta <- read.csv('GEx-metadata-with-other-env-layers-v2.csv')
Herb<- read.csv('Meta_SEV2019_v2_with_body_size.csv')

All2<-All %>% 
  left_join(Meta)%>%
  left_join(Climate)%>%
  left_join(Herb) %>% 
  mutate(ALLC3=(C3+C3.)) %>% 
  mutate(ALLC4=(C4+C4.)) %>% 
  mutate(PhotoMix=abs(ALLC3-ALLC4)) 

write.csv(All2, 'community_difference_allmetrics_siteavg_09June2020.csv', row.names=F)
