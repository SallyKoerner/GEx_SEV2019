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
photopath <- read.csv('percent_photosynthetic_pathway.csv')
phylorealms <- read.csv('PhyloRealms.csv')%>%select(-X)
domRankChange <- read.csv('sitelevel_domspecies_rankchange.csv')%>%select(-X)

###merge together
compDiffSite <- domDiff%>%
  left_join(compDiff)%>%
  left_join(RACdiff)%>%
  left_join(domIDdiff)%>%
  left_join(CmaxDiff)%>%
  left_join(codomDiff)%>%
  left_join(photopath)%>%
  left_join(phylorealms)%>%
  left_join(domRankChange)%>%
  mutate(site_dom=(GDom+UDom)/2) #calculate average dominance between grazed and ungrazed areas

# write.csv(compDiffSite, 'community_difference_allmetrics_siteavg_04262019.csv', row.names=F)